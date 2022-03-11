#include <fstream>
#include <json/json.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/boundary_facets.h>

#include "../Stitch.h"
#include "ElasticIO.h"

bool loadElastic(const std::string& path, ElasticSetup& setup, ElasticState& state)
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
	if (!parse(setup.density, json["density"]))
	{
		std::cout << "missing density." << std::endl;
		return false;
	}
	if (!parse(setup.penaltyK, json["collision_penalty"]))
	{
		std::cout << "missing collision penalty." << std::endl;
		setup.penaltyK = 0;
	}
	if (!parse(setup.innerEta, json["collision_eta"]))
	{
		std::cout << "missing collision eta." << std::endl;
		setup.innerEta = 0.5 * setup.thickness;
	}
	if (!parse(setup.perturb, json["perturb"]))
	{
		setup.perturb = 0;
	}
	if (!parse(setup.pressure, json["pressure"]))
	{
		std::cout << "missing pressure coef." << std::endl;
		setup.pressure = 0;
	}
	if (!parse(setup.restFlat, json["rest_flat"]))
	{
		std::cout << "missing rest flat flag, set default to true" << std::endl;
		setup.restFlat = true;
	}

	if (!parse(setup.isNoeHookean, json["isNoeHookean"]))
	{
		std::cout << "missing noehookean flag, set default to false (stvk material)" << std::endl;
		setup.isNoeHookean = false;
	}
	if (!parse(setup.tensionField, json["isTFT"]))
	{
		std::cout << "missing TFT flag, set default to false (full elastic model)." << std::endl;
		setup.tensionField = false;
	}

	if (!parse(setup.bendingType, json["bending_type"]))
	{
		std::cout << "missing bending type, set default to elastic." << std::endl;
		setup.bendingType = "elastic";
	}

	if (json["gravity"].isNull())
	{
		setup.gravity.setZero();
	}
	else
	{
		for (int i = 0; i < 3; i++)
			parse(setup.gravity(i), json["gravity"][i]);
	}
	if (!parse(setup.maxStepSize, json["max_stepsize"]))
	{
		setup.maxStepSize = 1.0;
	}
	if (!parse(setup.framefreq, json["frame_frequency"]))
	{
		setup.framefreq = 1;
	}
	if (!parse(setup.numInterp, json["num_interpolation"]))
	{
		setup.numInterp = 1;
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

	Eigen::MatrixXi F;
	if (!parse(setup.initMeshPath, json["init_mesh"]))
	{
		std::cout << "missing the initial guess path, set rest state as the initial state" << std::endl;
		state.initialGuess = setup.restV;
		F = setup.restF;
	}
	else if (!igl::readOBJ(filePathPrefix + setup.initMeshPath, state.initialGuess, F))
	{
		std::cout << "invalid the initial guess path, set rest state as the initial state" << std::endl;
		state.initialGuess = setup.restV;
		F = setup.restF;
	}
	state.mesh = MeshConnectivity(F);
	setup.sff->initializeExtraDOFs(state.initialEdgeDOFs, state.mesh, state.initialGuess);
	
	if (!parse(setup.curMeshPath, json["cur_mesh"]))
	{
		std::cout << "missing the initial guess path, set current state as the initial state" << std::endl;
		state.curPos = state.initialGuess;
		state.curEdgeDOFs = state.initialEdgeDOFs;
	}
	else if (!igl::readOBJ(filePathPrefix + setup.curMeshPath, state.curPos, F))
	{
		std::cout << "invalid the initial guess path, set rest state as the initial state" << std::endl;
		state.curPos = state.initialGuess;
		state.curEdgeDOFs = state.initialEdgeDOFs;
	}
	else if ((F - state.mesh.faces()).norm() != 0)
	{
		std::cout << "current mesh has the different mesh connectivity with the initial guess" << std::endl;
		return false;
	}
	setup.computeVertArea(state.curPos.rows(), state.mesh.faces());

	std::string prefix;
	index = setup.initMeshPath.rfind(".");
	prefix = setup.initMeshPath.substr(0, index);

	Eigen::MatrixXd unstitchedV;
	Eigen::MatrixXi unstitchedF;
	std::string meshName = filePathPrefix + prefix + std::string("_unstitched.obj");
	bool isUnstitchedGuess = false;
	if (igl::readOBJ(meshName, unstitchedV, unstitchedF))
	{
		isUnstitchedGuess = true;
		std::cout << "having unstitched shape" << std::endl;
		if (unstitchedV.rows() != setup.restV.rows())
		{
			std::cout << "miss matching between unstitched 3D shape and 2D patterns" << std::endl;
			return false;
		}
	}

	std::cout << "build laplacian" << std::endl;
	Eigen::VectorXi newIndex;
	std::vector<Eigen::Vector3i> bnd_edges;

	if (!isUnstitchedGuess)	
	{
		Eigen::MatrixXd refV;
		Eigen::MatrixXi refF;
		int nverts;

		if (state.initialGuess.rows() < setup.restV.rows())// no unstitched mesh provided and we do stitch patterns, we just build laplacian based on current stitched mesh TODO: fixed this by getting the unstitched pattern.
		{
			std::cout << "Missing unstitched 3D guess, using the current guess shape for construction of laplacian." << std::endl;
			refV = state.initialGuess;
			refF = state.mesh.faces();
			nverts = state.initialGuess.rows();
		}
		else
		{
			refV = setup.restV;
			refF = setup.restF;
			nverts = setup.restV.rows();
		}

		Eigen::MatrixXi bndV;
		Eigen::MatrixXi bndF;
		Eigen::MatrixXi OppBndV;

		igl::boundary_facets(refF, bndV, bndF, OppBndV);
		bnd_edges.clear();
		for (int i = 0; i < bndV.rows(); i++)
		{
			bnd_edges.push_back(Eigen::Vector3i(bndV(i, 0), bndV(i, 1), refF(bndF(i), OppBndV(i))));
		}

		newIndex.resize(refV.rows());
		for (int i = 0; i < newIndex.rows(); i++)
		{
			newIndex(i) = i;
		}

		setup.computeLaplacian(refV, refF, bnd_edges, newIndex, nverts);
	}
	else
	{
		Eigen::MatrixXd V;
		stitchMeshes(unstitchedV, unstitchedF, V, F, bnd_edges, newIndex);
		if ((state.mesh.faces() - F).norm() != 0)	// TODO: a better stitching code is neccessary 
		{
			std::cout << "missing match happened between loaded guess and stithing shape, the vertices indices may changed after stitching." << std::endl;
			igl::writeOBJ(prefix + std::string("_guess_stitched_trouble.obj"), V, F);
			/*std::cout << "guessed shape faces : " << state.mesh.nFaces() << ", stitched shape faces: " << F.rows() << std::endl;
			return false;*/
			std::cout << "using stored stitched mesh" << std::endl;
			V = state.curPos;
		}
		setup.computeLaplacian(setup.restV, setup.restF, bnd_edges, newIndex, V.rows());
	}

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

	int nverts = state.curPos.rows();
	//clamped vertices
	setup.clampedDOFs.clear();

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
		if (!ss || vid < 0 || vid >= nverts)
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
				setup.clampedDOFs[3 * vid + j] = state.curPos(vid, j);
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

	state.curEdgeDOFs = state.initialEdgeDOFs;
	if (parse(setup.curEdgeDOFsPath, json["curedge_DOFs"]) && state.initialEdgeDOFs.size())
	{
		std::string edgeDOFsPath = filePathPrefix + setup.curEdgeDOFsPath;
		std::cout << "load edge dofs in: " << edgeDOFsPath << std::endl;
		std::ifstream efs(edgeDOFsPath);
		state.curEdgeDOFs = state.initialEdgeDOFs;
		if (!efs)
		{
			std::cout << "failed to load edge dofs file in: " << setup.curEdgeDOFsPath << std::endl;
		}
		else
		{
			for (int i = 0; i < state.curEdgeDOFs.size(); i++)
			{
				efs >> state.curEdgeDOFs(i);
			}
		}
		std::cout << "edge dofs norm = " << state.curEdgeDOFs.norm() << std::endl;
	}
	return true;
}

bool saveElastic(const std::string& path, const ElasticSetup& setup, const ElasticState& state)
{
	Json::Value json;
	json["poisson_ratio"] = setup.PoissonsRatio;
	json["thickness"] = setup.thickness;
	json["youngs_modulus"] = setup.YoungsModulus;
	json["density"] = setup.density;
	json["collision_penalty"] = setup.penaltyK;
	json["pressure"] = setup.pressure;
	json["collision_eta"] = setup.innerEta;
	json["perturb"] = setup.perturb;
	json["rest_flat"] = setup.restFlat;
	json["isNeoHookean"] = setup.isNoeHookean;
	json["isTFT"] = setup.tensionField;
	json["bending_type"] = setup.bendingType;
	json["max_stepsize"] = setup.maxStepSize;
	json["gravity"].resize(3);
	for (int i = 0; i < 3; i++)
		json["gravity"][i] = setup.gravity(i);

	json["frame_frequency"] = setup.framefreq;
	json["num_interpolation"] = setup.numInterp;

	json["sff_type"] = setup.sffType;

	json["rest_mesh"] = setup.restMeshPath;
	json["init_mesh"] = setup.initMeshPath;
	json["cur_mesh"] = setup.curMeshPath;
	json["obs_mesh"] = setup.obstaclePath;

	json["curedge_DOFs"] = setup.curEdgeDOFsPath;
	json["clamped_DOFs"] = setup.clampedDOFsPath;

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

	Json::StyledStreamWriter writer;
	writer.write(file, json);
	file.close();

	if (!igl::writeOBJ(filePathPrefix + setup.curMeshPath, state.curPos, state.mesh.faces()))
		return false;

	if (state.curEdgeDOFs.size() > 0)
	{
		std::ofstream efs(filePathPrefix + setup.curEdgeDOFsPath);
		if (!efs)
			return false;
		else
		{
			for (int i = 0; i < state.curEdgeDOFs.size(); i++)
			{
				efs << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << state.curEdgeDOFs(i) << std::endl;
			}
		}
	}
	return true;
}