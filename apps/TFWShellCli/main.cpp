#include <CLI/CLI.hpp>
#include <igl/writeOBJ.h>
#include <regex>

#include "../../ShellSolver.h"
#include "../../TFWShell/TFWIO.h"
#include "../../ElasticShell/ElasticIO.h"

using namespace TFW;
TFWSetup setup;
TFWState curState;
TFWOptimizationParams TFWOptParams;
bool reInitialization = false;
bool quietOpt = false;
int upsampledTimes = 3;

std::string workingFolder;
std::string inputPath;
std::string outputFolder = "";

bool loadProblem(std::string loadPath)
{
    int index = loadPath.rfind("/");
	workingFolder = loadPath.substr(0, index);
	std::cout << "working folder: " << workingFolder << std::endl;

    bool ok = loadTFW(loadPath, setup, curState);
    return ok;
}

bool saveProblem(std::string saveFolder)
{
	mkdir(saveFolder);
	std::string wrinkledMeshName = std::regex_replace(setup.restMeshPath, std::regex(".obj"), "_wrinkledMesh.obj");
	igl::writeOBJ(saveFolder + "/" + wrinkledMeshName, curState.wrinkledPos, curState.wrinkledF);

	bool ok = saveTFW(saveFolder + "/data.json", setup, curState);
    return ok;
}

int main(int argc, char* argv[])
{
    CLI::App app("Tension Field + Wrinkles Simulator");
    app.add_option("input,-i,--input", inputPath, "Input model (json file)")->required()->check(CLI::ExistingFile);
	app.add_option("-o,--output", outputFolder, "Output Folder");
	app.add_option("-n,--numIter", TFWOptParams.iterations, "Number of iterations, default is 1000");
	app.add_option("-g,--gradTol", TFWOptParams.gradNorm, "The tolerance for gradient norm termination, default is 1e-6");
	app.add_option("-x,--xTol", TFWOptParams.xDelta, "The tolerance of variable update termination, default is 0");
	app.add_option("-f,--fTol", TFWOptParams.fDelta, "The tolerance of function update termination, default is 0");
	app.add_option("-u,--upsampleTimes", upsampledTimes, "The upsampling (subdivision) times to extract wrinkles, default is 3");
	app.add_flag("-q,--quiet", quietOpt, "Do Not print the optimization log, default is false");
	app.add_flag("-r,--reinitialization", reInitialization, "Reinitialize current amplitude and frequency based on the strain, default is false");

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
		outputFolder = workingFolder;  // use the input folder as the output one
	}

	// optimize for amplitude and frequency
	TFWOptParams.printLog = !quietOpt;
	if(reInitialization)
	{
		curState.reinitializeWrinkleVaribles(setup);
	}
	ShellSolver::TFWSQPSolver(setup, curState, workingFolder + std::regex_replace(setup.restMeshPath, std::regex(".obj"), ""), TFWOptParams);
	curState.getWrinkleMesh(setup, upsampledTimes);

	saveProblem(outputFolder);

    return 0;

}
