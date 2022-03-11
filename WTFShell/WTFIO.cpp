#include <fstream>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>

#include "WTFIO.h"
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

bool WTF::loadWTF(const std::string& path, WTFSetup& setup, WTFState& state)
{
	Json::Value json;
	Json::Reader reader;

	std::string filePath = path;
	std::replace(filePath.begin(), filePath.end(), '\\', '/'); // handle the backslash issue for windows

	int index = filePath.rfind("/");
	std::string filePathPrefix = filePath.substr(0, index + 1);

	std::ifstream file(path);

	bool parsingSuccessful = reader.parse(file, json);
	if (!parsingSuccessful) {
		std::cout << "Error reading file: " << path << std::endl;
		return false;
	}
	file.close();
	if (!parse(setup.PoissonsRatio, json["poisson_ratio"]))
	{
		std::cout << "missing poisson ratio." << std::endl;
		return false;
	}
	if (!parse(setup.thickness, json["thickness"]))
	{
		std::cout << "missing thickness." << std::endl;
		return false;
	}
	if (!parse(setup.YoungsModulus, json["youngs_modulus"]))
	{
		std::cout << "missing youngs modulus." << std::endl;
		return false;
	}

	if (!parse(setup.nasoqEps, json["nasoq_eps"]))
	{
		std::cout << "set default nasoq eps to 1e-6." << std::endl;
		setup.nasoqEps = 1e-6;
	}

	if (!parse(setup.quadNum, json["quadpoints"]))
	{
		std::cout << "missing quadrature points number, set default to 3." << std::endl;
		setup.quadNum = 3;
	}
	
	setup.buildQuadraturePoints();

	if (!parse(setup.restFlat, json["rest_flat"]))
	{
		std::cout << "missing the flag for whether the rest shape is flat, set default to true." << std::endl;
		setup.restFlat = true;
	}

	if (!parse(setup.clampedChosenVerts, json["clamp"]))
	{
		std::cout << "missing the flag for whether to clampe the given vertices, set default to true." << std::endl;
		setup.clampedChosenVerts = true;
	}
	
	if (!parse(setup.sffType, json["sff_type"]))
	{
		std::cout << "missing the sff type setting, set default to midedgetan." << std::endl;
		setup.sffType = "midedgeTan";
	}

	if (std::strcmp(setup.sffType.c_str(), "midedgeTan") == 0)
		setup.sff = std::make_shared<MidedgeAngleTanFormulation>();
	else if (std::strcmp(setup.sffType.c_str(), "midedgeSin") == 0)
		setup.sff = std::make_shared<MidedgeAngleSinFormulation>();
	else if (std::strcmp(setup.sffType.c_str(), "midedgeAve") == 0)
		setup.sff = std::make_shared<MidedgeAverageFormulation>();
	else
		setup.sff = std::make_shared<MidedgeAngleTanFormulation>();

	if (!parse(setup.restMeshPath, json["rest_mesh"]))
	{
		std::cout << "missing the rest mesh path." << std::endl;
		return false;
	}
	if (!igl::readOBJ(filePathPrefix + setup.restMeshPath, setup.restV, setup.restF))
		return false;
	setup.buildRestFundamentalForms();

	if (!parse(setup.obstaclePath, json["obs_mesh"]))
	{
		std::cout << "no obstacles" << std::endl;
	}
	Eigen::MatrixXd Vobs;
	Eigen::MatrixXi Fobs;
	if (igl::readOBJ(filePathPrefix + setup.obstaclePath, Vobs, Fobs))
	{
		Obstacle obstacle(Vobs, Fobs);
		setup.obs.push_back(obstacle);
	}

	if (!parse(setup.baseMeshPath, json["base_mesh"]))
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

	// TODO: modify this
	if (!parse(setup.clampedDOFsPath, json["clamped_DOFs"]))
	{
		std::cout << "missing clamped DOFs path." << std::endl;
	}
	std::ifstream ifs(filePathPrefix + setup.clampedDOFsPath);
	if (!ifs)
	{
		std::cout << "Missing " << setup.clampedDOFsPath << std::endl;
		return false;
	}
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

	int nverts = state.basePos.rows();
	int nedges = state.baseMesh.nEdges();

	std::string phiPath, dphiPath, ampPath;
	bool isAmpLoaded = parse(ampPath, json["amp_path"]);
	bool isDphiLoaded = parse(dphiPath, json["dphi_path"]);
	bool isPhiLoaded = parse(phiPath, json["phi_path"]);

	isAmpLoaded = isAmpLoaded ? parseTXT(filePathPrefix + ampPath, nverts, state.amplitude) : false;
	isDphiLoaded = isDphiLoaded ? parseTXT(filePathPrefix + dphiPath, nedges, state.dphi) : false;
	isPhiLoaded = isPhiLoaded ? parseTXT(filePathPrefix + phiPath, nverts, state.phi) : false;


	if (!isAmpLoaded || !isDphiLoaded || !isPhiLoaded)
	{
		std::cout << "amp file load status (0 for fail, 1 for succeed): " << isAmpLoaded << std::endl;
		std::cout << "phi file load status (0 for fail, 1 for succeed): " << isPhiLoaded << std::endl;
		std::cout << "dphi file load status (0 for fail, 1 for succeed): " << isDphiLoaded << std::endl;

		state.reinitializeWrinkleVaribles(setup);
	}
	return true;
}
bool WTF::saveWTF(const std::string& path, const WTFSetup& setup, WTFState& state)
{
	Json::Value json;
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

	std::ofstream file(path.c_str());
	if (!file)
	{
		std::cout << "output path doesn't exist." << std::endl;
		return false;
	}

	Json::StyledStreamWriter writer;
	writer.write(file, json);
	file.close();
	
	std::string filePath = path;
	std::replace(filePath.begin(), filePath.end(), '\\', '/'); // handle the backslash issue for windows

	int index = filePath.rfind("/");
	std::string filePathPrefix = filePath.substr(0, index + 1);

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