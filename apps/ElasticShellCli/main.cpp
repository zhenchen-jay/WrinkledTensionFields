#include <igl/hausdorff.h>
#include <iostream>
#include <random>
#include <CLI/CLI.hpp>
#include <regex>

#include "../../MeshLib/MeshConnectivity.h"
#include "../../ShellSolver.h"
#include "../../ElasticShell/ElasticIO.h"

int bendingType;

ElasticSetup setup;
ElasticState curState;
int numSteps;
FullSimOptimizationParams fullSimOptParams;

std::string workingFolder;
std::string inputPath;
std::string outputFolder = "";


SFFType sfftype = MidedgeAverage;
double perturbMag = 0;
bool quietOpt = false;


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

void setParameters()
{
    numSteps = setup.numInterp;
    if (setup.sffType == "midedgeAve")
        sfftype = MidedgeAverage;
    else if (setup.sffType == "midedgeSin")
        sfftype = MidedgeSin;
    else
        sfftype = MidedgeTan;

}

bool loadProblem(std::string loadPath)
{
	int index = loadPath.rfind("/");
	workingFolder = loadPath.substr(0, index);
	std::cout << "working folder: " << workingFolder << std::endl;

    bool ok = loadElastic(loadPath, setup, curState);
    
    if(ok)
    {
        setParameters();
    }
        
    return ok;
}

bool saveProblem(std::string saveFolder)
{
    mkdir(saveFolder);
    std::cout << "make dir done" << std::endl;
    bool ok = saveElastic(saveFolder + "/data.json", setup, curState);
    return ok;
}

int main(int argc, char* argv[])
{
    CLI::App app("Quasi-static Simulator");
    app.add_option("input,-i,--input", inputPath, "Input model (json file)")->required()->check(CLI::ExistingFile);
	app.add_option("-o,--output", outputFolder, "Output folder");
	app.add_option("-n,--numIter", fullSimOptParams.iterations, "Number of iterations, default is 1000");
	app.add_option("-g,--gradTol", fullSimOptParams.gradNorm, "The tolerance for gradient norm termination, default is 1e-6");
	app.add_option("-x,--xTol", fullSimOptParams.xDelta, "The tolerance of variable update termination, default is 0");
	app.add_option("-f,--fTol", fullSimOptParams.fDelta, "The tolerance of function update termination, default is 0");
	app.add_flag("-q,--quiet", quietOpt, "Do Not print the optimization log, default is false");
	app.add_option("-p,--randomPerturb", perturbMag, "add some random perturbation, default is 0 (not add)");

    try
    {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError& e)
    {
        return app.exit(e);
    }

	loadProblem(inputPath);
	if(outputFolder == "")
	{
		outputFolder = workingFolder;
	}

	jitter(perturbMag);
	for (int k = 1; k <= numSteps; k++) // how many interpolation steps should take (for clamping vertices)
	{
		ShellSolver::fullSimNewtonStaticSolver(setup, curState, workingFolder + std::regex_replace(setup.restMeshPath, std::regex(".obj"), ""), fullSimOptParams);
	}

	saveProblem(outputFolder);

    return 0;

}
