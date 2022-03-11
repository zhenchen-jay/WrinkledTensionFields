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
#include <igl/boundary_loop.h>

#include <Eigen/CholmodSupport>
#include <imgui/imgui.h>
#include <json/json.h>

#include "../../ShellSolver.h"
#include "../../Obstacle.h"

#include "../../ElasticShell/ElasticEnergy.h"

#include "../../CommonFunctions.h"
#include "../../PressureEnergy.h"

#include "../../EigenNASOQ.h"
#include "../../external/eigengurobi/Gurobi.h"
#include "../../ElasticShell/ElasticIO.h"


static StretchingType stretchingType;
int bendingType;
static ElasticMaterialType materialType;

ElasticSetup setup;
ElasticState curState;
int numSteps;
FullSimOptimizationParams fullSimOptParams;

double thickness;
double pressure;
double PoissonRatio;
double YoungsModulus;
double materialDensity;
double x_externalForce;
double y_externalForce;
double z_externalForce;
int frameFreq;
double noiseNorm;
double penaltyK;
double maxStepSize;
double innerEta;
int iterations;

bool isVisualizeOriginalMesh; // visualize the origin mesh
bool isVisualizeObstacle; // visualize the obstacle
bool restFlat; // plot tension line
double scale;   // plot scale
bool isSPDProj;

std::string filePathPrefix;
std::string defaultPath;
std::string elasticPath;


void liftVerts(double mag) // lifing the verts in direction z
{
	for (int i = 0; i < curState.curPos.rows(); i++)
		curState.curPos(i, 2) += mag;
}

void jitter(double magnitude)
{
	assert(magnitude >= 0);
	if (magnitude == 0)
		return;
	std::cout << "Noise magnitude: " << magnitude << std::endl;
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0, magnitude);

	Eigen::MatrixXd VN;

	igl::per_vertex_normals(curState.curPos, curState.mesh.faces(), VN);

	for (int i = 0; i < curState.curPos.rows(); i++)
	{
		if (setup.clampedDOFs.find(3 * i) == setup.clampedDOFs.end())
		{
			double offsetNoise = distribution(generator);
			curState.curPos.row(i) += offsetNoise * VN.row(i);
		}
	}
}

void reset()
{
	std::cout << std::endl << "Reset" << std::endl << std::endl;
	//curState.resetState();

}

void repaint(igl::opengl::glfw::Viewer& viewer)
{
	viewer.data().clear();

	/*viewer.data().set_mesh(curState.curPos, curState.mesh.faces());
	Eigen::MatrixXd colors(curState.curPos.rows(), 3);*/

	Eigen::MatrixXd curV, obsV;
	Eigen::MatrixXd curColor, obsColor;
	Eigen::MatrixXi curF, obsF;

	curV.resize(0, 3);
	curF.resize(0, 3);
	curColor.resize(0, 3);

	obsV.resize(0, 3);
	obsF.resize(0, 3);
	obsColor.resize(0, 3);
	

	bool isFaceBased = true;
	if (isVisualizeOriginalMesh)
	{
		curV = curState.curPos;
		curF = curState.mesh.faces();
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
	
	if (isVisualizeObstacle)
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

	std::cout << "Setting data (V, F): " << std::endl;
	Eigen::MatrixXd dataV, dataColors;
	Eigen::MatrixXi dataF;

	int totalNumV = curV.rows() + obsV.rows();
	int totalNumF = curF.rows() + obsF.rows();
	int totalNumColors = curColor.rows() + obsColor.rows();

	dataV.resize(totalNumV, 3);
	dataF.resize(totalNumF, 3);
	dataColors.resize(totalNumColors, 3);

	dataV.block(0, 0, curV.rows(), 3) = curV;
	dataV.block(curV.rows(), 0, obsV.rows(), 3) = obsV;
	
	dataF.block(0, 0, curF.rows(), 3) = curF;
	
	for (int i = 0; i < obsF.rows(); i++)
	{
		for(int j = 0; j < 3; j++)
			dataF(curF.rows() + i, j) = obsF(i, j) + curV.rows();
	}


	std::cout << "Setting Colors: " << std::endl;
	dataColors.block(0, 0, curColor.rows(), 3) = curColor;
	dataColors.block(curColor.rows(), 0, obsColor.rows(), 3) = obsColor;
	
	// reset data
	std::cout << "Reset Viewer" << std::endl;
	viewer.data().clear();
	viewer.data().set_face_based(isFaceBased);
	
	viewer.data().set_mesh(dataV, dataF);
	viewer.data().set_colors(dataColors);

	std::cout << "Reset Viewer done" << std::endl;
}

void setParameters()
{
	thickness = setup.thickness;
	materialDensity = setup.density;
	YoungsModulus = setup.YoungsModulus;
	PoissonRatio = setup.PoissonsRatio;
	numSteps = setup.numInterp;
	frameFreq = setup.framefreq;
	x_externalForce = setup.gravity(0);
	y_externalForce = setup.gravity(1);
	z_externalForce = setup.gravity(2);
	scale = 1;
	stretchingType = setup.tensionField ? tensionField : elasticStretching;
	iterations = 2000;
	//set to false for default instead
	bendingType = setup.bendingType == "elastic" ? elasticBending : noBending;
	if (setup.bendingType == "quadratic")
		bendingType = quadraticBending;
	noiseNorm = 0.0;
	penaltyK = setup.penaltyK;

	isSPDProj = setup.penaltyK > 0 ? false : true;

	pressure = setup.pressure;
	restFlat = setup.restFlat;
	innerEta = setup.innerEta;
	maxStepSize = setup.maxStepSize;

	isVisualizeOriginalMesh = true;
	isVisualizeObstacle = false;
}

bool sequentialRotateCylinder(const std::string& path)
{
    std::cout << "path: " << path << std::endl;
    int index = path.rfind("/");
    filePathPrefix = path.substr(0, index);
    std::cout << "filePathPrefix: " << filePathPrefix << std::endl;
    bool ok = loadElastic(path, setup, curState);
    if(!ok)
    {
        std::cout << "error in the loading cylinder file!" << std::endl;
        return ok;
    }

    setup.tensionField = true;
    setup.bendingType = BendingType::quadraticBending;  // use quad bending

    double zmax = curState.initialGuess.col(2).maxCoeff();
    double zmin = curState.initialGuess.col(2).minCoeff();

    std::vector<std::vector<int>> bnds;
    igl::boundary_loop(curState.mesh.faces(), bnds);

    int upbnds = 0;
    if(curState.initialGuess(bnds[0][0], 2) < (zmin + zmax) / 2)
        upbnds = 1;

    setParameters();

    curState.curPos = curState.initialGuess;    // set to the initial cylinder without twisted
    
    for(int i = 0; i < 20; i++)
    {
        double theta = (i + 1) * 2.0 / 180.0 * M_PI;
        std::cout << "rotation angle: " <<  (i + 1) * 2 << std::endl;
        Eigen::Matrix2d J;
        J << std::cos(theta), -std::sin(theta), std::sin(theta), std::cos(theta);

        for(auto& id: bnds[upbnds])
        {
            Eigen::Vector2d pos2d;
            pos2d<< curState.initialGuess(id, 0), curState.initialGuess(id, 1);
            pos2d = J * pos2d;
            setup.clampedDOFs[3 * id + 0] = pos2d(0);
            setup.clampedDOFs[3 * id + 1] = pos2d(1);
            setup.clampedDOFs[3 * id + 2] = curState.initialGuess(id, 2);
        }
        for(auto& id: bnds[1 - upbnds])
        {
            for(int j = 0; j < 3; j++)
                setup.clampedDOFs[3 * id + j] = curState.initialGuess(id, j);
        }

        for (int k = 1; k <= numSteps; k++)
        {
            ShellSolver::fullSimNewtonStaticSolver(setup, curState, filePathPrefix, fullSimOptParams);
        }

        igl::writeOBJ(filePathPrefix + "/rotation/" + std::to_string((i + 1) * 2) + "/cylinder_simulated.obj", curState.curPos, curState.mesh.faces());
        // write the clamped information
        std::string clampedFile = filePathPrefix + "/rotation/" + std::to_string((i + 1) * 2) + "/cylinder_clamped_vertices.dat";
        std::ofstream cfs = std::ofstream (clampedFile);
        if(cfs)
        {
            cfs << bnds[0].size() + bnds[1].size() << "\n";
            cfs << "#" << "\n";
            for(auto& id: bnds[0])
            {
                cfs << id << " " << curState.curPos.row(id) << "\n";
            }
            for(auto& id: bnds[1])
            {
                cfs << id << " " << curState.curPos.row(id) << "\n";
            }
        }

    }
    return true;
}

int main(int argc, char* argv[])
{
    sequentialRotateCylinder("/home/zchen96/Downloads/cylinder688verts/cylinder_Elastic.json");
	defaultPath = "../../Sims/disc_elastic/disc_elastic.json";
	if (argc >= 2)
		defaultPath = argv[1];
	std::cout << defaultPath << std::endl;
	int index = defaultPath.rfind("_");
	filePathPrefix = defaultPath.substr(0, index);
	// generateElasticJson(filePathPrefix);
	bool ok = loadElastic(defaultPath, setup, curState);
	int loadingTimes = 0;
	while (!ok)
	{
		defaultPath = "../" + defaultPath;
		ok = loadElastic(defaultPath, setup, curState);
		loadingTimes++;
		if (loadingTimes >= 5)
		{
			std::cout << "Failed to load Problem with the path" << defaultPath << std::endl;
			return -1;
		}
	}
	index = defaultPath.rfind("_");
	filePathPrefix = defaultPath.substr(0, index);
	elasticPath = defaultPath;
	std::cout << "Model Name is: " << filePathPrefix << std::endl;


	std::cout << "Set Parameters" << std::endl;
	// parameters loaded from the file
	setParameters();
	std::cout << "curEdgeDOFs norm = " << curState.curEdgeDOFs.norm() << std::endl;

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
					elasticPath = defaultPath;
				}
				else
				{
					int index = path.rfind("_");
					filePathPrefix = path.substr(0, index);
					std::cout << "Model Name is: " << filePathPrefix << std::endl;
					elasticPath = path;
				}
				ElasticSetup newSetup;
				ElasticState newState;
				bool ok = loadElastic(elasticPath, newSetup, newState);
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
					bool ok = saveElastic(path, setup, curState);
					if (!ok)
						std::cout << "Failed to save problem." << std::endl;
				}
			}
			if (ImGui::Button("Save Simed Mesh", ImVec2(-1, 0)))
			{
				const std::string path = igl::file_dialog_save();
				if (!igl::writeOBJ(path, curState.curPos, curState.mesh.faces()))
				{
					std::cout << "Failed to save simulated mesh." << std::endl;
				}
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
			//            ImGui::PopItemWidth();
		}
		// Draw options
		if (ImGui::CollapsingHeader("Draw Options", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (ImGui::InputDouble("plot scale", &scale))
			{
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
			float w = ImGui::GetContentRegionAvailWidth();
			float p = ImGui::GetStyle().FramePadding.x;
			if (ImGui::InputDouble("Thickness", &thickness))
			{
				if (thickness > 0)
				{
					setup.thickness = thickness;
				}
			}
			if (ImGui::InputDouble("Density", &materialDensity))
			{
				if (materialDensity > 0)
				{
					setup.density = materialDensity;
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
			if (ImGui::InputDouble("fext_x", &x_externalForce))
			{
				setup.gravity(0) = x_externalForce;	
			}
			if (ImGui::InputDouble("fext_y", &y_externalForce))
			{
				
				setup.gravity(1) = y_externalForce;
			}
			if (ImGui::InputDouble("fext_z", &z_externalForce))
			{
				setup.gravity(2) = z_externalForce;
			}

			if (ImGui::InputDouble("Collision Penalty", &penaltyK))
			{
					setup.penaltyK = penaltyK >= 0.0 ? penaltyK : 2000;
			}
			if (ImGui::InputDouble("Collision Eta", &innerEta))
			{
					setup.innerEta = innerEta >= 0.0 ? innerEta : 0.5 * thickness;
			}

			if (ImGui::InputDouble("Pressure", &pressure))
			{
				if (pressure > 0)
				{
					setup.pressure = pressure;
				}
			}
			if (ImGui::InputDouble("Noise Norm", &noiseNorm))
			{
				if (noiseNorm < 0)
				{
					noiseNorm = thickness / 2.0;
				}
			}
			if (ImGui::InputDouble("Max Step Size", &maxStepSize))
			{
				if (maxStepSize > 0)
				{
					setup.maxStepSize = maxStepSize;
				}
			}
		}
		if (ImGui::CollapsingHeader("Optimization Parameters", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (ImGui::InputInt("Frame Frequence", &frameFreq))
			{
				if (frameFreq > 0)
				{
					setup.framefreq = frameFreq;
				}
			}
			if (ImGui::InputInt("Interpolation Steps", &numSteps))
			{
				if (numSteps > 0)
				{
					setup.numInterp = numSteps;
				}
				fullSimOptParams.interp = 1.0 / numSteps;
			}

			if (ImGui::InputInt("num of iterations", &iterations))
			{
				if (!(iterations > 0))
				{
					fullSimOptParams.iterations = 2000;
				}
				else
				{
					fullSimOptParams.iterations = iterations;
				}
			}
			if (ImGui::InputDouble("Grad Tol", &fullSimOptParams.gradNorm))
			{
				if (!(fullSimOptParams.gradNorm >= 0))
				{
					fullSimOptParams.gradNorm = 1e-6;
				}
			}
			if (ImGui::InputDouble("delta_f", &fullSimOptParams.fDelta))
			{
				if (!(fullSimOptParams.fDelta >= 0))
				{
					fullSimOptParams.fDelta = 1e-10;
				}
			}
			if (ImGui::InputDouble("delta_x", &fullSimOptParams.xDelta))
			{
				if (!(fullSimOptParams.xDelta >= 0))
				{
					fullSimOptParams.xDelta = 1e-10;
				}
			}


			ImGui::Checkbox("SPD Proj H", &fullSimOptParams.isProjH);
			ImGui::Checkbox("StVK TBB Parallel", &fullSimOptParams.isParallel);
		}

		if (ImGui::Combo("Elastic Material", (int*)&materialType, "StVK\0NoeHookean\0\0"))
		{
			if (materialType == NeoHookean)
				setup.isNoeHookean = true;
			else
				setup.isNoeHookean = false;
		}
		if (ImGui::CollapsingHeader("Elastic Energy Type", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (ImGui::Combo("Stretching", (int*)&stretchingType, "Elastic\0TensionField\0\0"))
			{
				switch (stretchingType)
				{
				case elasticStretching:
					setup.tensionField = false;
					break;
				case tensionField:
					setup.tensionField = true;
					break;
				default:
					break;
				}
			}
			if (ImGui::Combo("Bending", (int*)&bendingType, "Elastic\0Quadratic\0NoBending\0\0"))
			{
				switch (bendingType)
				{
				case elasticBending:
					setup.bendingType = "elastic";
					break;
				case quadraticBending:
					setup.bendingType = "quadratic";
					break;
				case noBending:
					setup.bendingType = "nobending";
					break;
				default:
					break;
				}
			}
		}
		if (ImGui::Button("Optimize ElasticShell", ImVec2(-1, 0)))
		{
			double reg = 1e-6;
			int funcEvals = 0;
			double interp = 0;
			double prevEnergy = 0;
			double curEnergy = 0;

			/*double _lameAlpha = setup.YoungsModulus * setup.PoissonsRatio / (1.0 - setup.PoissonsRatio * setup.PoissonsRatio);
			double _lameBeta = setup.YoungsModulus / 2.0 / (1.0 + setup.PoissonsRatio);
			std::shared_ptr<ElasticShellMaterial> mat;
			mat = std::make_shared<StVKMaterial>();

			std::cout << "stretching = " << elasticStretchingEnergy(curState.mesh, curState.curPos, curState.curEdgeDOFs, _lameAlpha, _lameBeta, setup.thickness, setup.abars, *sff, *mat, NULL, NULL, false, true) << std::endl;
			std::cout << "bending = " << elasticBendingEnergy(curState.mesh, curState.curPos, curState.curEdgeDOFs, _lameAlpha, _lameBeta, setup.thickness, setup.abars, setup.bbars, *sff, *mat, NULL, NULL, false, true) << std::endl;*/
			jitter(noiseNorm);
			for (int k = 1; k <= numSteps; k++)
			{	
				ShellSolver::fullSimNewtonStaticSolver(setup, curState, filePathPrefix, fullSimOptParams);
			}
			saveElastic(elasticPath, setup, curState);
			repaint(viewer);
		}

		if (ImGui::Button("Reset ElasticShell", ImVec2(-1, 0)))
		{
			reset();
			curState.resetState();
			repaint(viewer);
		}

		if (ImGui::CollapsingHeader("Useful functions"))
		{
			if (ImGui::Button("Compute Hausdorff Distance", ImVec2(-1, 0)))
			{
				Eigen::MatrixXd V1, V2;
				Eigen::MatrixXi F1, F2;
				std::string fileName1 = igl::file_dialog_open();
				if (!igl::readOBJ(fileName1, V1, F1))
				{
					std::cout << "Failed to load first mesh" << std::endl;
				}
				else
				{
					std::string fileName2 = igl::file_dialog_open();
					if (!igl::readOBJ(fileName2, V2, F2))
					{
						std::cout << "Failed to load second mesh" << std::endl;
					}
					else
					{
						double d = 0;
						igl::hausdorff(V1, F1, V2, F2, d);
						std::cout << "Hausdorff Distance: " << d << std::endl;
					}
				}
			}

			if (ImGui::Button("convert json file", ImVec2(-1, 0)))
			{
				const std::string path = igl::file_dialog_save();
				int index = path.rfind(".");
				std::string myprefix = path.substr(0, index);
				generateElasticJson(myprefix);
			}
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
