#include <igl/hausdorff.h>
#include <iostream>
#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>

#include <CLI/CLI.hpp>

#include "polyscope/polyscope.h"
#include "polyscope/messages.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/view.h"
#include "polyscope/point_cloud.h"


#include "../../ShellSolver.h"
#include "../../TFWShell/TFWIO.h"

using namespace TFW;
TFWSetup setup;
TFWState curState;
TFWState initialState;

int quadType;
TFWOptimizationParams TFWOptParams;

double thickness;
double pressure;
double PoissonRatio;
double YoungsModulus;

bool isVisualizeWrinkles; // visualize interior points along the wrinkles
bool isVisualizeObstacle; // visualize the obstacle
int upsamplingTimes = 2;
double scale;   // plot scale
bool isUseV1Term;
bool isUseV2Term;
bool isFixBnd;
RoundingType roundingType;


std::string filePathPrefix;
std::string defaultPath;
std::string TFWPath;

float surfShiftx = 0;
float surfShifty = 0;
float surfShiftz = 0;

float preShiftx = 0;
float preShifty = 0;
float preShiftz = 0;

float boundx = 0;
float boundy = 0;
float boundz = 0;


void reset()
{
	std::cout << std::endl << "Reset" << std::endl << std::endl;
	curState = initialState;
}


void setParameters()
{
	thickness = setup.thickness;
	YoungsModulus = setup.YoungsModulus;
	PoissonRatio = setup.PoissonsRatio;
	scale = 1;

	isVisualizeWrinkles = false;
	isVisualizeObstacle = false;
	quadType = 1;
	upsamplingTimes = 4;
	initialState = curState;

	isUseV1Term = false;
	isUseV2Term = true;

	roundingType = CWFRound;
	isFixBnd = false;
}

bool loadProblem(std::string loadPath = "")
{
    if(loadPath == "")
        defaultPath = igl::file_dialog_open();
    else
        defaultPath = loadPath;
    
    int index = defaultPath.rfind("_");
    filePathPrefix = defaultPath.substr(0, index);
    TFWPath = defaultPath;
    std::cout << "Model Name is: " << filePathPrefix << std::endl;
    bool ok = loadTFW(defaultPath, setup, curState);
    
    if(ok)
    {
        setParameters();
        boundx = curState.basePos.col(0).maxCoeff() - curState.basePos.col(0).minCoeff();
        boundy = curState.basePos.col(1).maxCoeff() - curState.basePos.col(1).minCoeff();
        boundz = curState.basePos.col(2).maxCoeff() - curState.basePos.col(2).minCoeff();
    }
        
    return ok;
}

bool saveProblem(std::string savePath = "")
{
    if (savePath == "")
    {
        savePath = igl::file_dialog_save();
    }
	if(!isVisualizeWrinkles)
	{
		std::cout << "construct wrinkle mesh" << std::endl;
		curState.getWrinkleMesh(setup, upsamplingTimes, isUseV1Term, isUseV2Term, roundingType, isFixBnd);
	}

    bool ok = saveTFW(savePath, setup, curState);
    return ok;
}

void updateView()
{
    polyscope::registerSurfaceMesh("base mesh", curState.basePos, curState.baseMesh.faces());
    polyscope::getSurfaceMesh("base mesh")->setSurfaceColor({ 80 / 255.0, 122 / 255.0, 91 / 255.0 });

	if(isVisualizeWrinkles)
	{
		curState.getWrinkleMesh(setup, upsamplingTimes, isUseV1Term, isUseV2Term, roundingType, isFixBnd);
		polyscope::registerSurfaceMesh("wrinkled mesh", curState.wrinkledPos, curState.wrinkledF);
		polyscope::getSurfaceMesh("wrinkled mesh")->setSurfaceColor({ 122 / 255.0, 80 / 255.0, 91 / 255.0 });
	}
	else
	{
		if(polyscope::hasSurfaceMesh("wrinkled mesh"))
			polyscope::getSurfaceMesh("wrinkled mesh")->setEnabled(false);
	}

	if(isVisualizeObstacle)
	{
		if(!setup.obs.empty())
		{
			for(int i = 0; i < setup.obs.size(); i++)
			{
				auto initObs = polyscope::registerSurfaceMesh("init obstacle " + std::to_string(i), setup.obs[i].V, setup.obs[i].F);
				auto curObjs = polyscope::registerSurfaceMesh("current obstacle " + std::to_string(i), setup.obs[i].V, setup.obs[i].F);

				initObs->setSurfaceColor({1, 1, 1});
				curObjs->setSurfaceColor({1, 1, 1});
			}
		}
	}
	else
	{
		if(!setup.obs.empty())
		{
			for(int i = 0; i < setup.obs.size(); i++)
			{
				if(polyscope::hasSurfaceMesh("init obstacle " + std::to_string(i)))
				{
					auto initObs = polyscope::getSurfaceMesh("init obstacle " + std::to_string(i));
					initObs->setEnabled(false);
				}
				if(polyscope::hasSurfaceMesh("current obstacle " + std::to_string(i)))
				{
					auto curObjs = polyscope::getSurfaceMesh("current obstacle " + std::to_string(i));
					curObjs->setEnabled(false);
				}

			}
		}
	}
}


void callback() {
    ImGui::PushItemWidth(100);
    float w = ImGui::GetContentRegionAvailWidth();
    float p = ImGui::GetStyle().FramePadding.x;
    if (ImGui::Button("Load", ImVec2((w - p) / 2.f, 0)))
    {
        loadProblem();
		updateView();
    }
    ImGui::SameLine(0, p);
    if (ImGui::Button("Save", ImVec2((w - p) / 2.f, 0)))
    {
        saveProblem();
    }

    // Draw options
    if (ImGui::CollapsingHeader("Draw Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::Checkbox("Visualize Obstacle", &isVisualizeObstacle))
        {
            updateView();
        }
	    if (ImGui::Checkbox("Visualize Wrinkles", &isVisualizeWrinkles))
	    {
			updateView();
		}
    }

	// upsampling options
	if (ImGui::CollapsingHeader("Upsampling Options", ImGuiTreeNodeFlags_DefaultOpen))
	{
		if (ImGui::InputInt("Upsampling Times", &upsamplingTimes))
		{
			updateView();
		}
		if (ImGui::Checkbox("Use V1 Term", &isUseV1Term))
		{
			updateView();
		}
		if (ImGui::Checkbox("Use V2 Term", &isUseV2Term))
		{
			updateView();
		}
		if (ImGui::Checkbox("Fix bnd", &isFixBnd))
		{
			updateView();
		}
		if (ImGui::Combo("Round Type", (int*)&roundingType, "Comiso\0CWFRound\0\0"))
		{
			updateView();
		}

	}

    if (ImGui::CollapsingHeader("wrinkled mesh shift", ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::DragFloat("x", &surfShiftx, 0.05, 0, 4))
        {
            if(polyscope::hasSurfaceMesh("wrinkled mesh"))
            {
                polyscope::getSurfaceMesh("wrinkled mesh")->translate({ (surfShiftx - preShiftx) * boundx, 0, 0});
                if(!setup.obs.empty())
                {
                    for(int i = 0; i < setup.obs.size(); i++)
                    {
	                    if(polyscope::hasSurfaceMesh("current obstacle " + std::to_string(i)))
	                    {
		                    polyscope::getSurfaceMesh("current obstacle " + std::to_string(i))->translate({ (surfShiftx - preShiftx) * boundx, 0, 0});
						}

                    }
                }
                preShiftx = surfShiftx;
            }
        }

        if (ImGui::DragFloat("y", &surfShifty, 0.05, 0, 4))
        {
            if(polyscope::hasSurfaceMesh("wrinkled mesh"))
            {
                polyscope::getSurfaceMesh("wrinkled mesh")->translate({ 0, (surfShifty - preShifty) * boundy, 0 });
                if(!setup.obs.empty())
                {
                    for(int i = 0; i < setup.obs.size(); i++)
                    {
	                    if(polyscope::hasSurfaceMesh("current obstacle " + std::to_string(i)))
	                    {
		                    polyscope::getSurfaceMesh("current obstacle " + std::to_string(i))->translate({ 0, (surfShifty - preShifty) * boundy, 0 });
						}

                    }
                }
                preShifty = surfShifty;
            }
        }

        if (ImGui::DragFloat("z", &surfShiftz, 0.05, 0, 4))
        {
            if(polyscope::hasSurfaceMesh("wrinkled mesh"))
            {
                polyscope::getSurfaceMesh("wrinkled mesh")->translate({ 0, 0, (surfShiftz - preShiftz) * boundz });
                if(!setup.obs.empty())
                {
                    for(int i = 0; i < setup.obs.size(); i++)
                    {
						if(polyscope::hasSurfaceMesh("current obstacle " + std::to_string(i)))
						{
							polyscope::getSurfaceMesh("current obstacle " + std::to_string(i))->translate({ 0, 0, (surfShiftz - preShiftz) * boundz });
						}

                    }
                }
                preShiftz = surfShiftz;
            }
        }
    }

    if (ImGui::CollapsingHeader("Material Parameters", ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::InputDouble("Thickness", &thickness))
        {
            if (thickness > 0)
            {
                setup.thickness = thickness;
            }
        }
        if (ImGui::InputDouble("Young's Modulus", &YoungsModulus))
        {
            if (YoungsModulus > 0)
            {
                setup.YoungsModulus = YoungsModulus;
            }
        }
        if (ImGui::InputDouble("Poisson's Ratio", &PoissonRatio))
        {
            if (PoissonRatio >= 0)
            {
                setup.PoissonsRatio = PoissonRatio;
            }
        }
    }
    if (ImGui::CollapsingHeader("Optimization Parameters", ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::InputInt("num of iterations", &TFWOptParams.iterations))
        {
            if (!(TFWOptParams.iterations > 0))
            {
	            TFWOptParams.iterations = 2000;
            }
        }
        if (ImGui::InputDouble("Grad Tol", &TFWOptParams.gradNorm))
        {
            if (!(TFWOptParams.gradNorm >= 0))
            {
	            TFWOptParams.gradNorm = 1e-6;
            }
        }
        if (ImGui::InputDouble("delta_f", &TFWOptParams.fDelta))
        {
            if (!(TFWOptParams.fDelta >= 0))
            {
	            TFWOptParams.fDelta = 1e-10;
            }
        }
        if (ImGui::InputDouble("delta_x", &TFWOptParams.xDelta))
        {
            if (!(TFWOptParams.xDelta >= 0))
            {
	            TFWOptParams.xDelta = 1e-10;
            }
        }
        ImGui::Checkbox("StVK TBB Parallel", &TFWOptParams.isParallel);
    }

    if (ImGui::CollapsingHeader("TFW Type", ImGuiTreeNodeFlags_DefaultOpen))
    {
	    if (ImGui::Combo("Quadrature number", &quadType, "one\0three\0six\0sixteen\0\0"))
	    {
		    switch (quadType)
		    {
			    case 0:
				    setup.quadNum = 1;
				    break;
			    case 1:
				    setup.quadNum = 3;
				    break;
			    case 2:
				    setup.quadNum = 6;
				    break;
			    case 3:
				    setup.quadNum = 16;
				    break;
			    default:
				    setup.quadNum = 3;
				    break;
		    }
		    setup.buildQuadraturePoints();
		    std::cout << "number of quadrature points used : " << setup.quadPoints.size() << std::endl;
	    }

	    if(ImGui::Checkbox("clamp chosen verts", &setup.clampedChosenVerts))
        {
            std::vector<Eigen::Vector3d> clampedPts;

	        std::unordered_set<int> clampedVids;

	        for (auto &it : setup.clampedDOFs)
	        {
		        clampedVids.insert(it.first / 3);
	        }
			for(auto & it: clampedVids)
			{
				clampedPts.push_back(curState.basePos.row(it).segment<3>(0).transpose());
			}
	        auto pscloud = polyscope::registerPointCloud("clamped pts", clampedPts);
        }

	    if (ImGui::Button("Reinitalization", ImVec2(-1, 0)))
	    {
		    curState.reinitializeWrinkleVaribles(setup);
		    initialState = curState;
		    updateView();
	    }

        if (ImGui::Button("Optimize TFW", ImVec2(-1, 0)))
        {
			ShellSolver::TFWSQPSolver(setup, curState, filePathPrefix, TFWOptParams);
            updateView();
        }
    }   
    ImGui::PopItemWidth();
}

int main(int argc, char* argv[])
{
    CLI::App app("Tension Field + Wrinkles Simulator");
    app.add_option("input,-i,--input", defaultPath, "Input model (json file)")->check(CLI::ExistingFile);

    try
    {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError& e)
    {
        return app.exit(e);
    }
    
    // Options
//	polyscope::options::autocenterStructures = true;
    polyscope::view::windowWidth = 1024;
    polyscope::view::windowHeight = 1024;

    // Initialize polyscope
    polyscope::init();

    polyscope::view::upDir = polyscope::view::UpDir::ZUp;

    // Add the callback
    polyscope::state::userCallback = callback;

    polyscope::options::groundPlaneHeightFactor = 0.25; // adjust the plane height

    if(loadProblem(defaultPath))
        updateView();

    // Show the gui
    polyscope::show();

    return 0;

}
