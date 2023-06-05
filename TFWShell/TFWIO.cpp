#include <fstream>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <json.hpp>

#include "TFWIO.h"
#include "PhiEstimate.h"

bool parseTXT(const std::string& path, const int& n, Eigen::VectorXd& v)
{
	std::ifstream vfile(path);
	if (!vfile)
		return false;
	v.resize(n);
	v.setZero();
	for (int i = 0; i < n; i++)
		vfile >> v(i);
	return true;
}

bool TFW::loadTFW(const std::string& path, TFWSetup& setup, TFWState& state)
{
	std::string filePath = path;
	std::replace(filePath.begin(), filePath.end(), '\\', '/'); // handle the backslash issue for windows

	int index = filePath.rfind("/");
	std::string filePathPrefix = filePath.substr(0, index + 1);

	std::ifstream inputJson(filePath);

	if (!inputJson)
	{
		std::cerr << "json file does not exist: " << filePath << std::endl;
		return false;
	}

	nlohmann::json jval;
	inputJson >> jval;

	setup.PoissonsRatio = jval["poisson_ratio"];
	setup.thickness = jval["thickness"];
	setup.YoungsModulus = jval["youngs_modulus"];

	setup.nasoqEps = 1e-6;
	if(jval.contains(std::string_view{"nasoq_eps"}))
	{
		setup.nasoqEps = jval["nasoq_eps"];
	}

	setup.quadNum = 3;
	if(jval.contains(std::string_view{"quadpoints"}))
	{
		setup.quadNum = jval["quadpoints"];
	}
	setup.buildQuadraturePoints();

	if (jval.contains(std::string_view{ "rest_flat" }))
	{
		setup.restFlat = jval["rest_flat"];
	}
	else
		setup.restFlat = true;

	if (jval.contains(std::string_view{ "clamp" }))
	{
		setup.clampedChosenVerts = jval["clamp"];
	}
	else
		setup.clampedChosenVerts = true;

	if (jval.contains(std::string_view{ "sff_type" }))
	{
		setup.sffType = jval["sff_type"];
	}
	else
		setup.sffType = "midedgeTan";

	if (std::strcmp(setup.sffType.c_str(), "midedgeTan") == 0)
		setup.sff = std::make_shared<MidedgeAngleTanFormulation>();
	else if (std::strcmp(setup.sffType.c_str(), "midedgeSin") == 0)
		setup.sff = std::make_shared<MidedgeAngleSinFormulation>();
	else if (std::strcmp(setup.sffType.c_str(), "midedgeAve") == 0)
		setup.sff = std::make_shared<MidedgeAverageFormulation>();
	else
		setup.sff = std::make_shared<MidedgeAngleTanFormulation>();

	if (jval.contains(std::string_view{ "rest_mesh" }))
	{
		setup.restMeshPath = jval["rest_mesh"];
		if (!igl::readOBJ(filePathPrefix + setup.restMeshPath, setup.restV, setup.restF))
			return false;
	}
	else
	{
		std::cout << "missing the rest mesh path." << std::endl;
		return false;
	}
	setup.buildRestFundamentalForms();

	if (jval.contains(std::string_view{ "obs_mesh" }))
	{
		setup.obstaclePath = jval["obs_mesh"];
		Eigen::MatrixXd Vobs;
		Eigen::MatrixXi Fobs;
		if (igl::readOBJ(filePathPrefix + setup.obstaclePath, Vobs, Fobs))
		{
			Obstacle obstacle(Vobs, Fobs);
			setup.obs.push_back(obstacle);
		}
	}

	if (jval.contains(std::string_view{ "base_mesh" }))
	{
		setup.baseMeshPath = jval["base_mesh"];
	}
	else
	{
		std::cout << "missing the base mesh path." << std::endl;
		return false;
	}
	Eigen::MatrixXi baseF;
	if (!igl::readOBJ(filePathPrefix + setup.baseMeshPath, state.basePos, baseF))
		return false;
	state.baseMesh = MeshConnectivity(baseF);
	setup.sff->initializeExtraDOFs(state.baseEdgeDOFs, state.baseMesh, state.basePos);
	locatePotentialPureTensionFaces(setup.abars, state.basePos, state.baseMesh.faces(), state.tensionFaces);

	std::cout << "compute the curvature info of the base mesh." << std::endl;
	state.computeBaseCurvature();           // This step is important. 

	if(!jval.contains(std::string_view{ "clamped_DOFs" }))
	{
		std::cout << "missing clamped DOFs file." << std::endl;
	}
	else
	{
		setup.clampedDOFsPath = jval["clamped_DOFs"];
		std::ifstream ifs(filePathPrefix + setup.clampedDOFsPath);
		if (!ifs)
		{
			std::cout << "Missing " << setup.clampedDOFsPath << std::endl;
		}
		else
		{
			int nclamped;
			ifs >> nclamped;
			char dummy;
			ifs >> dummy;
			if (!ifs)
			{
				std::cout << "Error in " << setup.clampedDOFsPath << std::endl;
				return false;
			}
			ifs.ignore(std::numeric_limits<int>::max(), '\n');
			std::cout << "num of clamped DOFs: " << nclamped << std::endl;
			for (int i = 0; i < nclamped; i++)
			{
				std::string line;
				std::getline(ifs, line);
				std::stringstream ss(line);

				int vid;
				ss >> vid;
				if (!ss || vid < 0 || vid >= state.basePos.rows())
				{
					std::cout << "Error in " << setup.clampedDOFsPath << std::endl;
					return false;
				}
				std::string x; // x, y, z
				ss >> x;
				if (!ss)
				{
					std::cout << "-using rest position for clamped vertices " << vid << std::endl;
					for (int j = 0; j < 3; j++)
					{
						setup.clampedDOFs[3 * vid + j] = state.basePos(vid, j);
					}
				}
				else
				{
					for (int j = 0; j < 3; j++)
					{
						if (x[0] != '#')
						{
							setup.clampedDOFs[3 * vid + j] = std::stod(x);
						}
						ss >> x;
					}
				}
			}

			if (!ifs)
			{
				std::cout << "Error in " << setup.clampedDOFsPath << std::endl;
				return false;
			}
		}
	}

	int nverts = state.basePos.rows();
	int nedges = state.baseMesh.nEdges();

	int id = setup.restMeshPath.rfind(".");

	std::string modelName = setup.restMeshPath.substr(0, id);

	setup.ampPath = modelName + "_amplitude_simulated.txt";
	setup.dphiPath = modelName + "_dphi_simulated.txt";
	setup.phiPath = modelName + "_phi_simulated.txt";

	if(jval.contains(std::string_view{"amp_path"}))
	{
		setup.ampPath = jval["amp_path"];
	}
	if(jval.contains(std::string_view{"dphi_path"}))
	{
		setup.dphiPath = jval["dphi_path"];
	}
	if(jval.contains(std::string_view{"phi_path"}))
	{
		setup.phiPath = jval["phi_path"];
	}

	bool isAmpLoaded = parseTXT(filePathPrefix + setup.ampPath, nverts, state.amplitude);
	bool isDphiLoaded = parseTXT(filePathPrefix + setup.dphiPath, nedges, state.dphi);
	bool isPhiLoaded = parseTXT(filePathPrefix + setup.phiPath, nverts, state.phi);


	if (!isAmpLoaded || !isDphiLoaded || !isPhiLoaded)
	{
		std::cout << "amp file load status (0 for fail, 1 for succeed): " << isAmpLoaded << std::endl;
		std::cout << "phi file load status (0 for fail, 1 for succeed): " << isPhiLoaded << std::endl;
		std::cout << "dphi file load status (0 for fail, 1 for succeed): " << isDphiLoaded << std::endl;

		state.reinitializeWrinkleVaribles(setup);
	}
	return true;
}
bool TFW::saveTFW(const std::string& path, const TFWSetup& setup, TFWState& state)
{
	nlohmann::json json;
	json["poisson_ratio"] = setup.PoissonsRatio;
	json["thickness"] = setup.thickness;
	json["youngs_modulus"] = setup.YoungsModulus;
	json["nasoq_eps"] = setup.nasoqEps;
	json["quadpoints"] = setup.quadNum;
	json["rest_flat"] = setup.restFlat;
	json["clamp"] = setup.clampedChosenVerts;
	json["sff_type"] = setup.sffType;

	json["rest_mesh"] = setup.restMeshPath;
	json["base_mesh"] = setup.baseMeshPath;
	json["clamped_DOFs"] = setup.clampedDOFsPath;
	json["obs_mesh"] = setup.obstaclePath;

	json["dphi_path"] = setup.dphiPath;
	json["amp_path"] = setup.ampPath;
	json["phi_path"] = setup.phiPath;

	std::string filePath = path;
	std::replace(filePath.begin(), filePath.end(), '\\', '/'); // handle the backslash issue for windows

	int index = filePath.rfind("/");
	std::string filePathPrefix = filePath.substr(0, index + 1);

	std::ofstream file(path.c_str());
	if (!file)
	{
		std::cout << "output path doesn't exist." << path << std::endl;
		return false;
	}
	file << std::setw(4) << json << std::endl;
	std::cout << "save file in: " << path << std::endl;

	// phi 
	if (state.phi.size())
	{
		std::string phiFileName = filePathPrefix + setup.phiPath;
		std::ofstream pfs(phiFileName);
		if (!pfs)
			return false;
		for (int i = 0; i < state.phi.size(); i++)
		{
			pfs << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << state.phi(i) << std::endl;
		}

	}

	// dphi 
	if (state.dphi.size())
	{
		std::string dphiFileName = filePathPrefix + setup.dphiPath;
		std::ofstream dpfs(dphiFileName);
		if (!dpfs)
			return false;
		for (int i = 0; i < state.dphi.size(); i++)
		{
			dpfs << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << state.dphi(i) << std::endl;
		}
	}


	// amplitude
	if (state.amplitude.size())
	{
		std::string amplitudeFileName = filePathPrefix + setup.ampPath;
		std::ofstream afs(amplitudeFileName);
		if (!afs)
			return false;
		afs.precision(15);
		for (int i = 0; i < state.amplitude.size(); i++)
		{
			afs << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << state.amplitude(i) << std::endl;
		}
	}

	return true;
}