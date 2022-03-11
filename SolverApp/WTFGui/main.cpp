#include <igl/opengl/glfw/Viewer.h>
#include <igl/writePLY.h>
#include <igl/hausdorff.h>
#include "../../MeshConnectivity.h"
#include <iostream>
#include <random>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/signed_distance.h>
#include <igl/cotmatrix.h>
#include <igl/per_vertex_normals.h>

#include <Eigen/CholmodSupport>
#include <imgui/imgui.h>
#include <json/json.h>

#include "../../SecondFundamentalForm/SecondFundamentalFormDiscretization.h"
#include "../../SecondFundamentalForm/MidedgeAngleSinFormulation.h"
#include "../../SecondFundamentalForm/MidedgeAngleTanFormulation.h"
#include "../../SecondFundamentalForm/MidedgeAverageFormulation.h"
#include "../../ShellSolver.h"
#include "../../WTFShell/MeshUpsampling.h"
#include "../../Obstacle.h"
#include "../../WTFShell/PhiEstimate.h"

#include "../../WTFShell/WTFShell.h"

#include "../../CommonFunctions.h"

#include "../../ZuenkoFormula.h"
#include "../../WTFShell/WTFIO.h"


WTFSetup setup;
WTFState curState;
WTFState initialState;

int quadType;
WTFOptimizationParams WTFOptParams;

double thickness;
double pressure;
double PoissonRatio;
double YoungsModulus;

bool isVisualizeWrinkles; // visualize interior points along the wrinkles
bool isVisualizeObstacle; // visualize the obstacle
int upsamplingTimes = 2;
bool isVisualizeOriginalMesh;
bool clampedChosenVerts;
bool restFlat; // plot tension line
double scale;   // plot scale
bool isSPDProj;
bool isShowProblemFaces; // tension faces after smoothing
bool isUseV1Term;
bool isUseV2Term;


std::string filePathPrefix;
std::string defaultPath;
std::string WTFPath;



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

	isSPDProj = true;
	quadType = 1;	// 
	
	isVisualizeOriginalMesh = true;
	upsamplingTimes = 4;
	initialState = curState;


	//clampedChosenVerts = !setup.clampedChosenVerts ? setup.clampedChosenVerts : true;
	clampedChosenVerts = setup.clampedChosenVerts;

	isShowProblemFaces = false;
	isUseV1Term = false;
	isUseV2Term = true;
}

void repaint(igl::opengl::glfw::Viewer& viewer)
{
	viewer.data().clear();
	Eigen::MatrixXd curV, obsV, wrinkledV;
	Eigen::MatrixXd curColor, obsColor, wrinkledColor;
	Eigen::MatrixXi curF, obsF, wrinkledF;

	// edges
	Eigen::MatrixXd PStart, PEnd;
	Eigen::RowVector3d edgeColor(0, 0, 0);

	curV.resize(0, 3);
	curF.resize(0, 3);
	curColor.resize(0, 3);

	obsV.resize(0, 3);
	obsF.resize(0, 3);
	obsColor.resize(0, 3);

	wrinkledV.resize(0, 3);
	wrinkledF.resize(0, 3);
	wrinkledColor.resize(0, 3);

	PStart.resize(0, 3);
	PEnd.resize(0, 3);
	

	bool isFaceBased = true;
	if (isVisualizeOriginalMesh)
	{
		curV = curState.basePos;
		curF = curState.baseMesh.faces();
		std::cout << "visualize base mesh." << std::endl;
		std::cout << "num of verts: " << curV.rows() << ", num of faces: " << curF.rows() << std::endl;

		isFaceBased = false;
		curColor.resize(curV.rows(), 3);
		curColor.col(0).setConstant(1.0);
		curColor.col(1).setConstant(1.0);
		curColor.col(2).setConstant(0.0);
		for (auto& it : setup.clampedDOFs)
		{
			int vid = it.first / 3;
			curColor(vid, 0) = 1.0;
			curColor(vid, 1) = 0.0;
			curColor(vid, 2) = 0.0;
		}
	}
	
	if (isVisualizeObstacle && setup.obs.size())
	{
		std::cout << "visualize obstacle: " << std::endl;
		obsV = setup.obs[0].V;
		obsF = setup.obs[0].F;

		int colorRows = isFaceBased ? obsF.rows() : obsV.rows();
		obsColor.resize(colorRows, 3);
		for (int i = 0; i < colorRows; i++)
		{
			obsColor.row(i) << 1, 1, 1;
		}
	}
	
	if (isVisualizeWrinkles)
	{
		curState.getWrinkleMesh(setup, upsamplingTimes, isUseV1Term, isUseV2Term);
		igl::writeOBJ(filePathPrefix + "_wrinkledMesh.obj", curState.wrinkledPos, curState.wrinkledF);

		wrinkledV = curState.wrinkledPos;
		wrinkledF = curState.wrinkledF;
		int colorRows = isFaceBased ? wrinkledF.rows() : wrinkledV.rows();
		wrinkledColor.resize(colorRows, 3);
		for (int i = 0; i < colorRows; i++)
		{
			wrinkledColor.row(i) << 1, 0, 0;
		}
	}

	Eigen::MatrixXd cutStartP, cutEndP;
	Eigen::Vector3d cutcolor;
	if (isShowProblemFaces)
	{
		Eigen::Vector3d vertcolor;
		vertcolor << 1, 1, 1;
		std::set<int> problemFaces = curState.tensionFaces;

		isFaceBased = true;

		curColor.resize(curState.baseMesh.nFaces(), 3);
		curColor.col(0).setConstant(1.0);
		curColor.col(1).setConstant(1.0);
		curColor.col(2).setConstant(0.0);

		for (auto& fid : problemFaces)
		{
			curColor.row(fid) << vertcolor;
		}

		std::cout << "num of problem faces: " << problemFaces.size() << std::endl;
	}

	std::cout << "Setting data (V, F): " << std::endl;
	Eigen::MatrixXd dataV, dataColors;
	Eigen::MatrixXi dataF;

	int totalNumV = curV.rows() + obsV.rows() + wrinkledV.rows();
	int totalNumF = curF.rows() + obsF.rows() + wrinkledF.rows();
	int totalNumColors = curColor.rows() + obsColor.rows() + wrinkledColor.rows();

	dataV.resize(totalNumV, 3);
	dataF.resize(totalNumF, 3);
	dataColors.resize(totalNumColors, 3);

	dataV.block(0, 0, curV.rows(), 3) = curV;
	dataV.block(curV.rows(), 0, obsV.rows(), 3) = obsV;
	dataV.block(curV.rows() + obsV.rows(), 0, wrinkledV.rows(), 3) = wrinkledV;

	dataF.block(0, 0, curF.rows(), 3) = curF;
	
	for (int i = 0; i < obsF.rows(); i++)
	{
		for(int j = 0; j < 3; j++)
			dataF(curF.rows() + i, j) = obsF(i, j) + curV.rows();
	}

	for (int i = 0; i < wrinkledF.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
			dataF(curF.rows() + obsF.rows() + i, j) = wrinkledF(i, j) + curV.rows() + obsV.rows();
	}

	std::cout << "Setting Colors: " << std::endl;
	dataColors.block(0, 0, curColor.rows(), 3) = curColor;
	dataColors.block(curColor.rows(), 0, obsColor.rows(), 3) = obsColor;
	dataColors.block(curColor.rows() + obsColor.rows(), 0, wrinkledColor.rows(), 3) = wrinkledColor;

	// reset data
	std::cout << "Reset Viewer" << std::endl;
	viewer.data().clear();
	viewer.data().set_face_based(isFaceBased);
	
	viewer.data().set_mesh(dataV, dataF);
	viewer.data().set_colors(dataColors);

	std::cout << "adding edges: " << std::endl;
	if(PStart.rows())
		viewer.data().add_edges(PStart, PEnd, edgeColor);
	if(cutStartP.rows())
		viewer.data().add_edges(cutStartP, cutEndP, cutcolor);

	std::cout << "Reset Viewer done" << std::endl;
}

bool sequentialRotateCylinder(std::string path)
{
    WTFState lastState;
    for(int i = 0; i < 20; i++)
    {
        std::string mypath = path + "/rotation/" + std::to_string((i + 1) * 2) + "/cylinder_WTF.json";
        bool ok = loadWTF(mypath, setup, curState);
        if(!ok)
        {
            std::cout << "failed to load data in: " << mypath << std::endl;
            return false;
        }
        WTFPath = mypath;
        int index = mypath.rfind("_");
        filePathPrefix = mypath.substr(0, index);
        std::cout << "Model Name is: " << filePathPrefix << std::endl;

        std::cout << "Set Parameters" << std::endl;
        setParameters();

        curState.reinitializeWrinkleVaribles(setup);

        if(i != 0)
        {
            curState.phi = lastState.phi;
            curState.dphi = lastState.dphi;
            curState.amplitude = lastState.amplitude;
        }

        std::cout << "Start with: " << std::endl;
        std::cout << "norm of phi : " << curState.phi.norm() << std::endl;
        std::cout << "norm of dphi : " << curState.dphi.norm() << std::endl;
        std::cout << "norm of amplitude : " << curState.amplitude.norm() << std::endl;

        std::cout << "file path prefix: " << filePathPrefix << std::endl;

        ShellSolver::WTFSQPSolver(setup, curState, filePathPrefix, WTFOptParams);
        saveWTF(WTFPath, setup, curState);

        std::cout << "end with: " << std::endl;
        std::cout << "norm of phi : " << curState.phi.lpNorm<Eigen::Infinity>() << std::endl;
        std::cout << "norm of dphi : " << curState.dphi.lpNorm<Eigen::Infinity>() << std::endl;
        std::cout << "norm of sqrt amplitude : " << curState.amplitude.lpNorm<Eigen::Infinity>() << std::endl;

		std::string ampFile = path + "/rotation/amp/amp_" + std::to_string(i) + ".txt";
		std::ofstream afs(ampFile);
		afs << curState.amplitude << std::endl; 

		std::string dphiFile = path + "/rotation/omega/omega_" + std::to_string(i) + ".txt";
		std::ofstream dfs(dphiFile);
		dfs << curState.dphi << std::endl; 

		std::string meshFile = path + "/rotation/mesh/mesh_" + std::to_string(i) + ".obj";
		igl::writeOBJ(meshFile, curState.basePos, curState.baseMesh.faces());

        Eigen::MatrixXd NV, NVSeam;
        Eigen::MatrixXi NF, NFSeam;


        Eigen::MatrixXd curV = curState.basePos;
        Eigen::MatrixXi curF = curState.baseMesh.faces();
        Eigen::VectorXd curAmp = curState.amplitude;
        Eigen::VectorXd curPhi = curState.phi;


        std::string filename = filePathPrefix + "_wrinkledMesh.obj";
        igl::writeOBJ(filename.c_str(), curState.wrinkledPos, curState.wrinkledF);

         lastState = curState;
    }

    return true;
}

int main(int argc, char* argv[])
{
    sequentialRotateCylinder("/home/zchen96/Downloads/cylinder688verts/");
	defaultPath = "../../Sims/disc/disc_WTF.json";
	if (argc >= 2)
		defaultPath = argv[1];
	std::cout << defaultPath << std::endl;
	int index = defaultPath.rfind("_");
	filePathPrefix = defaultPath.substr(0, index);
	std::cout << "Model Name is: " << filePathPrefix << std::endl;
	// generateWTFJson(filePathPrefix);

	bool ok = loadWTF(defaultPath, setup, curState);
	int loadingTimes = 0;
	while (!ok)
	{
		defaultPath = "../" + defaultPath;
		ok = loadWTF(defaultPath, setup, curState);
		loadingTimes++;
		if (loadingTimes >= 5)
		{
			std::cout << "Failed to load Problem with the path" << defaultPath << std::endl;
			return -1;
		}
	}
	WTFPath = defaultPath;
	index = defaultPath.rfind("_");
	filePathPrefix = defaultPath.substr(0, index);
	std::cout << "Model Name is: " << filePathPrefix << std::endl;
	
	std::cout << "Set Parameters" << std::endl;
	setParameters();

	igl::opengl::glfw::Viewer viewer;

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);



	// Add content to the default menu window
	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		//        menu.draw_viewer_menu();
		if (ImGui::CollapsingHeader("I/O", ImGuiTreeNodeFlags_DefaultOpen))
		{
			float w = ImGui::GetContentRegionAvailWidth();
			float p = ImGui::GetStyle().FramePadding.x;
			if (ImGui::Button("Load Problem", ImVec2(ImVec2((w - p) / 2.f, 0))))
			{
				const std::string path = igl::file_dialog_open();
				if (path == "")
				{
					std::cout << "Failed to load problem. Set to the defualt path: " << defaultPath << std::endl;
					WTFPath = defaultPath;
				}
				else
				{
					int index = path.rfind("_");
					filePathPrefix = path.substr(0, index);
					std::cout << "Model Name is: " << filePathPrefix << std::endl;
					WTFPath = path;
				}
				WTFSetup newSetup;
				WTFState newState;
				bool ok = loadWTF(WTFPath, newSetup, newState);
				if (!ok)
				{
					std::cout << "Failed to load files !!!!" << std::endl;
					return;
				}
				else
				{
					setup = newSetup;
					curState = newState;
					setParameters();
					repaint(viewer);
					viewer.core().align_camera_center(viewer.data().V, viewer.data().F);
				}
			}
			ImGui::SameLine(0, p);
			if (ImGui::Button("Save Problem", ImVec2((w - p) / 2.f, 0)))
			{
				const std::string path = igl::file_dialog_save();
				if (path == "")
				{
					std::cout << "Invalid saving path" << std::endl;
				}
				else
				{
					bool ok = saveWTF(path, setup, curState);
					if (!ok)
						std::cout << "Failed to save problem." << std::endl;
				}
			}
			
			if (ImGui::Button("Save Wrinkled Mesh", ImVec2(-1, 0)))
			{

				std::string path = igl::file_dialog_save();
				Eigen::MatrixXd NV, NVSeam;
				Eigen::MatrixXi NF, NFSeam;
				Eigen::MatrixXd upsampledTFTV, soupPhiV, soupProblemV;
				Eigen::MatrixXi upsampledTFTF, soupPhiF, soupProblemF;
				Eigen::VectorXd upsampledAmp, soupPhi;


				Eigen::MatrixXd curV = curState.basePos;
				Eigen::MatrixXi curF = curState.baseMesh.faces();
				Eigen::VectorXd curAmp = curState.amplitude;
				Eigen::VectorXd curPhi = curState.phi;
				std::set<int> clampedVerts;
				if (setup.clampedChosenVerts)
				{
					for (auto& it : setup.clampedDOFs)
					{
						int vid = it.first / 3;
						if (clampedVerts.count(vid) == 0)
							clampedVerts.insert(vid);
					}
				}
				std::set<int> problemFaces;
				roundPhiFromDphiCutbyTension(curState.basePos, curState.baseMesh.faces(), curV, curF, setup.abars, curState.amplitude, curState.dphi, GurobiRound, curState.phi, curPhi, curAmp, problemFaces);
				wrinkledMeshUpsamplingUncut(curState.basePos, curState.baseMesh.faces(), setup.restV, setup.restF, curV, curF, problemFaces, clampedVerts,
					&NV, &NF, &upsampledTFTV, &upsampledTFTF, &soupPhiV, &soupPhiF, &soupProblemV, &soupProblemF, &upsampledAmp, &soupPhi,
					curAmp, curPhi, *setup.sff, setup.YoungsModulus, setup.PoissonsRatio, upsamplingTimes, Loop, isUseV1Term, isUseV2Term);

				// stitchMeshesWithTol(NVSeam, NFSeam, NV, NF, 1e-5);
				if (!igl::writeOBJ(path, NV, NF))
				{
					path = filePathPrefix + "_wrinkledMesh.obj";
					std::cout << "Invalid path, use current path: " << path << std::endl;
				}
				size_t lastindex = path.find_last_of(".");
				std::string rawname = path.substr(0, lastindex);
				std::string upsampledTFTpath = rawname + "_upsampledTFT.obj";
				igl::writeOBJ(upsampledTFTpath, upsampledTFTV, upsampledTFTF);
				std::string soupPath = rawname + "_phiSoup.obj";
				igl::writeOBJ(soupPath, soupPhiV, soupPhiF);
				std::string ampPath = rawname + "_upsampledAmp.csv";
				std::ofstream ampFile(ampPath);
				for (int i = 0; i < upsampledAmp.size(); i++)
					ampFile << upsampledAmp[i] << ",\t" << 3.14159 << std::endl; // fix for buggy Houdini import
				std::string phiPath = rawname + "_phiSoup.csv";
				std::ofstream phiFile(phiPath);
				for (int i = 0; i < soupPhi.size(); i++)
					phiFile << soupPhi[i] << ",\t" << 3.14159 << std::endl;

				std::string problempath = rawname + "_problem.obj";
				igl::writeOBJ(problempath, soupProblemV, soupProblemF);
			}
		}
		if (ImGui::CollapsingHeader("Viewing Options", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (ImGui::Button("Center object", ImVec2(-1, 0)))
			{
				viewer.core().align_camera_center(viewer.data().V, viewer.data().F);
			}
			if (ImGui::Button("Snap canonical view", ImVec2(-1, 0)))
			{
				viewer.snap_to_canonical_quaternion();
			}

			// Select rotation type
			int rotation_type = static_cast<int>(viewer.core().rotation_type);
			//int rotation_type = 0;
			static Eigen::Quaternionf trackball_angle = Eigen::Quaternionf::Identity();
			static bool orthographic = true;
			if (ImGui::Combo("Camera Type", &rotation_type, "Trackball\0Two Axes\02D Mode\0\0"))
			{
				using RT = igl::opengl::ViewerCore::RotationType;
				auto new_type = static_cast<RT>(rotation_type);
				if (new_type != viewer.core().rotation_type)
				{
					if (new_type == RT::ROTATION_TYPE_NO_ROTATION)
					{
						trackball_angle = viewer.core().trackball_angle;
						orthographic = viewer.core().orthographic;
						viewer.core().trackball_angle = Eigen::Quaternionf::Identity();
						viewer.core().orthographic = true;
					}
					else if (viewer.core().rotation_type == RT::ROTATION_TYPE_NO_ROTATION)
					{
						viewer.core().trackball_angle = trackball_angle;
						viewer.core().orthographic = orthographic;
					}
					viewer.core().set_rotation_type(new_type);
				}
			}

			// Orthographic view
			ImGui::Checkbox("Orthographic view", &(viewer.core().orthographic));
		}
		 // Draw options
		if (ImGui::CollapsingHeader("Draw Options", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (ImGui::InputDouble("plot scale", &scale))
			{
				repaint(viewer);
			}
			if (ImGui::Checkbox("Visualize Wrinkles", &isVisualizeWrinkles))
			{
				repaint(viewer);
			}

			if (ImGui::Checkbox("v1 in plane correction", &isUseV1Term))
			{
				repaint(viewer);
			}
			if (ImGui::Checkbox("v2 in plane correction", &isUseV2Term))
			{
				repaint(viewer);
			}
			if (ImGui::InputInt("upsampling times", &upsamplingTimes))
			{
				if (upsamplingTimes >= 0)
					repaint(viewer);
			}
			if (ImGui::Checkbox("Visualize Original Mesh", &isVisualizeOriginalMesh))
			{
				repaint(viewer);
			}
			if (ImGui::Checkbox("Visualize Obstacle", &isVisualizeObstacle))
			{
				repaint(viewer);
			}

			if (ImGui::Checkbox("show problem faces", &isShowProblemFaces))
			{
				repaint(viewer);
			}

			if (ImGui::Checkbox("Invert normals", &(viewer.data().invert_normals)))
			{
				viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_NORMAL;
			}
			ImGui::ColorEdit4("Background", viewer.core().background_color.data(),
				ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
			ImGui::ColorEdit4("Line color", viewer.data().line_color.data(),
				ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
			ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.3f);
			ImGui::DragFloat("Shininess", &(viewer.data().shininess), 0.05f, 0.0f, 100.0f);
			//ImGui::PopItemWidth();
		}

		// Overlays
		auto make_checkbox = [&](const char* label, unsigned int& option)
		{
			return ImGui::Checkbox(label,
				[&]() { return viewer.core().is_set(option); },
				[&](bool value) { return viewer.core().set(option, value); }
			);
		};

		if (ImGui::CollapsingHeader("Overlays", ImGuiTreeNodeFlags_DefaultOpen))
		{
			make_checkbox("Wireframe", viewer.data().show_lines);
			make_checkbox("Fill", viewer.data().show_faces);
			bool showvertid = viewer.data().show_vertex_labels != 0;
			if (ImGui::Checkbox("Show vertex labels", &showvertid))
			{
				viewer.data().show_vertex_labels = (showvertid ? 1 : 0);
			}
			bool showfaceid = viewer.data().show_face_labels != 0;
			if (ImGui::Checkbox("Show faces labels", &showfaceid))
			{
				viewer.data().show_face_labels = showfaceid;
			}
		}
	};
	menu.callback_draw_custom_window = [&]()
	{
		// Define next window position + size
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(0.0, 0.0), ImGuiCond_FirstUseEver);
		ImGui::Begin(
			"Optimization", nullptr,
			ImGuiWindowFlags_NoSavedSettings
		);

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
	
		if (ImGui::CollapsingHeader("Optimzation Parameters", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (ImGui::InputInt("Max iterations", &WTFOptParams.iterations))
			{
				if (!(WTFOptParams.iterations > 0))
				{
					WTFOptParams.iterations = 1e3;
				}
			}
			if (ImGui::InputDouble("Grad Tol", &WTFOptParams.gradNorm))
			{
				if (!(WTFOptParams.gradNorm >= 0))
				{
					WTFOptParams.gradNorm = 1e-10;
				}
			}
			if (ImGui::InputDouble("delta_f", &WTFOptParams.fDelta))
			{
				if (!(WTFOptParams.fDelta >= 0))
				{
					WTFOptParams.fDelta = 1e-10;
				}
			}
			if (ImGui::InputDouble("delta_x", &WTFOptParams.xDelta))
			{
				if (!(WTFOptParams.xDelta >= 0))
				{
					WTFOptParams.xDelta = 1e-10;
				}
			}
			ImGui::Checkbox("WTF TBB Parallel", &WTFOptParams.isParallel);
		}
		if (ImGui::CollapsingHeader("WTF Type", ImGuiTreeNodeFlags_DefaultOpen))
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
			
			ImGui::Checkbox("clamp chosen(red) verts", &setup.clampedChosenVerts);
			if (ImGui::Button("Reinitalization", ImVec2(-1, 0)))
			{
				curState.reinitializeWrinkleVaribles(setup);
				initialState = curState;
				repaint(viewer);
			}

		}

		if (ImGui::Button("Optimize WTF", ImVec2(-1, 0)))
		{
			std::cout << "Start with: " << std::endl;
			std::cout << "norm of phi : " << curState.phi.norm() << std::endl;
			std::cout << "norm of dphi : " << curState.dphi.norm() << std::endl;
			std::cout << "norm of amplitude : " << curState.amplitude.norm() << std::endl;

			std::cout << "file path prefix: " << filePathPrefix << std::endl;
			
			ShellSolver::WTFSQPSolver(setup, curState, filePathPrefix, WTFOptParams);
			saveWTF(WTFPath, setup, curState);

			std::cout << "end with: " << std::endl;
			std::cout << "norm of phi : " << curState.phi.lpNorm<Eigen::Infinity>() << std::endl;
			std::cout << "norm of dphi : " << curState.dphi.lpNorm<Eigen::Infinity>() << std::endl;
			std::cout << "norm of sqrt amplitude : " << curState.amplitude.lpNorm<Eigen::Infinity>() << std::endl;

			Eigen::MatrixXd NV, NVSeam;
			Eigen::MatrixXi NF, NFSeam;


			Eigen::MatrixXd curV = curState.basePos;
			Eigen::MatrixXi curF = curState.baseMesh.faces();
			Eigen::VectorXd curAmp = curState.amplitude;
			Eigen::VectorXd curPhi = curState.phi;

			
			std::string filename = filePathPrefix + "_wrinkledMesh.obj";
			igl::writeOBJ(filename.c_str(), curState.wrinkledPos, curState.wrinkledF);

			isVisualizeWrinkles = true;
			repaint(viewer);
		}
		ImGui::End();
	};
	viewer.data().set_face_based(false);
	std::cout << "set data done!" << std::endl;
	repaint(viewer);
	std::cout << "Repaint done!" << std::endl;
	viewer.launch();
	return 0;

}
