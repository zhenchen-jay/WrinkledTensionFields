#include <fstream>
#include <string>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/boundary_facets.h>
#include <json.hpp>

#include "ElasticIO.h"

bool loadElastic(const std::string& path, ElasticSetup& setup, ElasticState& state)
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
    setup.density = jval["density"];
    if (jval.contains(std::string_view{ "collision_penalty" }))
    {
        setup.penaltyK = jval["collision_penalty"];
    }
    else
        setup.penaltyK = 0;

    if (jval.contains(std::string_view{ "collision_eta" }))
    {
        setup.innerEta = jval["collision_eta"];
    }
    else
        setup.innerEta = 0.5 * setup.thickness;

    if (jval.contains(std::string_view{ "perturb" }))
    {
        setup.perturb = jval["perturb"];
    }
    else
        setup.perturb = 0;

    if (jval.contains(std::string_view{ "pressure" }))
    {
        setup.pressure = jval["pressure"];
    }
    else
        setup.pressure = 0;

    if (jval.contains(std::string_view{ "rest_flat" }))
    {
        setup.restFlat = jval["rest_flat"];
    }
    else
        setup.restFlat = true;

    if (jval.contains(std::string_view{ "stretching_type" }))
    {
        setup.strecthingType = jval["stretching_type"];
    }
    else
        setup.strecthingType = "StVK";

    if (jval.contains(std::string_view{ "bending_type" }))
    {
        setup.bendingType = jval["bending_type"];
    }
    else
        setup.bendingType = "elastic";

    if (jval.contains(std::string_view{ "gravity" }))
    {
        setup.gravity << jval["gravity"][0], jval["gravity"][1], jval["gravity"][2];
    }
    else
        setup.gravity.setZero();

    if (jval.contains(std::string_view{ "num_interpolation" }))
    {
        setup.numInterp = jval["num_interpolation"];
    }
    else
        setup.numInterp = 1;

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

    Eigen::MatrixXi F;
    if (!jval.contains(std::string_view{ "init_mesh" }))
    {
        std::cout << "missing the initial guess path, set rest state as the initial state" << std::endl;
        state.initialGuess = setup.restV;
        F = setup.restF;
    }
    else
    {
        setup.initMeshPath = jval["init_mesh"];
        if (!igl::readOBJ(filePathPrefix + setup.initMeshPath, state.initialGuess, F))
        {
            std::cout << "invalid the initial guess path, set rest state as the initial state" << std::endl;
            state.initialGuess = setup.restV;
            F = setup.restF;
        }
    }
    state.mesh = MeshConnectivity(F);
    setup.sff->initializeExtraDOFs(state.initialEdgeDOFs, state.mesh, state.initialGuess);

    if (!jval.contains(std::string_view{ "cur_mesh" }))
    {
        std::cout << "missing the current stat path, set current state as the initial state" << std::endl;
        state.curPos = state.initialGuess;
        state.curEdgeDOFs = state.initialEdgeDOFs;
    }
    else
    {
        setup.curMeshPath = jval["cur_mesh"];
        if (!igl::readOBJ(filePathPrefix + setup.curMeshPath, state.curPos, F))
        {
            std::cout << "invalid the current stat path, set rest state as the initial state" << std::endl;
            state.curPos = state.initialGuess;
            state.curEdgeDOFs = state.initialEdgeDOFs;
        }

    }
    if ((F - state.mesh.faces()).norm() != 0)
    {
        std::cout << "current mesh has the different mesh connectivity with the initial guess" << std::endl;
        return false;
    }
    setup.computeVertArea(state.curPos.rows(), state.mesh.faces());

    std::string prefix;
    index = setup.initMeshPath.rfind(".");
    prefix = setup.initMeshPath.substr(0, index);

    std::cout << "build laplacian" << std::endl;
    Eigen::VectorXi newIndex;
    std::vector<Eigen::Vector3i> bnd_edges;

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

    nverts = state.curPos.rows();
    //clamped vertices
    setup.clampedDOFs.clear();

    if (!jval.contains(std::string_view{ "clamped_DOFs" }))
    {
        std::cout << "missing clamped DOFs path." << std::endl;
    }
    else
    {
        setup.clampedDOFsPath = jval["clamped_DOFs"];
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
    }



    state.curEdgeDOFs = state.initialEdgeDOFs;
    if (jval.contains(std::string_view{ "curedge_DOFs" }) && state.initialEdgeDOFs.size())
    {
        setup.curEdgeDOFsPath = jval["curedge_DOFs"];
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
		if(setup.obs.size() && setup.penaltyK == 0)
			setup.penaltyK = 1e3;
	}

    return true;
}

bool saveElastic(const std::string& path, const ElasticSetup& setup, const ElasticState& state)
{
    nlohmann::json json;
    json["poisson_ratio"] = setup.PoissonsRatio;
    json["thickness"] = setup.thickness;
    json["youngs_modulus"] = setup.YoungsModulus;
    json["density"] = setup.density;
    json["collision_penalty"] = setup.penaltyK;
    json["pressure"] = setup.pressure;
    json["collision_eta"] = setup.innerEta;
    json["perturb"] = setup.perturb;
    json["rest_flat"] = setup.restFlat;
    json["stretching_type"] = setup.strecthingType;
    json["bending_type"] = setup.bendingType;
    json["max_stepsize"] = setup.maxStepSize;
    json["gravity"] = { setup.gravity(0), setup.gravity(1), setup.gravity(2) };
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
    file << std::setw(4) << json << std::endl;
    std::cout << "save file in: " << path << std::endl;

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