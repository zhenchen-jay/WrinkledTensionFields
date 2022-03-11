#include  <igl/boundary_loop.h>
#include <igl/writeOBJ.h>
#include <igl/readOBJ.h>
#include "ShellSolver.h"
#include "WTFShell/WTFShell.h"

void ShellSolver::WTFSQPSolver(const WTFSetup& setup, WTFState& curState, std::string filePrefix, const WTFOptimizationParams params)
{
	std::cout << "init amp norm: " << curState.amplitude.norm() << std::endl;
	std::cout << "init dphi norm: " << curState.dphi.norm() << std::endl;

	int nverts = curState.basePos.rows();
	int nedges = curState.baseMesh.nEdges();

	std::map<int, double> clampedDOFs;

	if (setup.clampedChosenVerts) // whether we clamped the vertices which were clamped in the TFT step
	{
		std::map<int, double>::const_iterator it;
		for (it = setup.clampedDOFs.begin(); it != setup.clampedDOFs.end(); it++)
		{
			curState.amplitude[it->first / 3] = 0;
			clampedDOFs[it->first / 3] = 0;
		}
	}

	if (curState.amplitude.minCoeff() < 0)
	{
		std::cout << "negative amplitude, use abs value instead!" << std::endl;
		curState.amplitude = curState.amplitude.cwiseAbs();
	}

	WTFModel model;
	model.initialization(setup, curState, &clampedDOFs, filePrefix, true, params.isParallel);

	Eigen::VectorXd initX;
	model.convertCurState2Variables(curState, initX);
	std::cout << initX.norm() << std::endl;

	auto C = model.buildIntegrabilityConstraints();
	//model.checkPD4ConstraintedHess(initX);
	double energy = model.value(initX);

	std::cout << "initial residual of integrabilty : " << (C * initX).norm() << std::endl;
	std::cout << "started energy: " << energy << std::endl;

	model.testGradientAndHessian(initX);
	model.testValueAndGradient(initX);

	cppoptlib::NASOQSolver<WTFModel> solver;
	cppoptlib::Criteria<double> s;
	s.iterations = params.iterations;
	s.fDelta = params.fDelta;
	s.gradNorm = params.gradNorm;
	s.xDelta = params.xDelta;
	solver.setStopCriteria(s);
	std::cout << "Stop Criteria: " << std::endl;
	std::cout << "Max iterations: " << s.iterations << ", ||g|| = " << s.gradNorm << ", ||delta_f|| = " << s.fDelta << ", ||delta_x|| = " << s.xDelta << std::endl;
	
	solver.minimize(model, initX);
	std::cout << "final residual of integrabilty : " << (C * initX).norm() << std::endl;

	model.convertVariables2CurState(initX, curState);
	curState.getWrinkleMesh(setup);
}

void ShellSolver::fullSimNewtonStaticSolver(const ElasticSetup& setup, ElasticState& curState, std::string filePrefix, const FullSimOptimizationParams params)
{
	ElasticShellModel model;
	bool ok = model.initialization(setup, curState, filePrefix, params.interp, params.isProjH, params.isParallel);
	if (!ok)
    {
        std::cout << "initialization failed." << std::endl;
        return;
    }

	Eigen::VectorXd initX;
	std::cout<<"initialization finished!"<<std::endl;
	Eigen::MatrixXi F;
	//igl::readOBJ("C:/Users/csyzz/Projects/tensionField/Sims/dress/MD/dress_simulated_1.obj", curState.curPos, F);
   
	model.convertCurState2Variables(curState, initX);
	std::cout<<"convertion finished!"<<std::endl;
	double energy = model.value(initX);

	Eigen::VectorXd grad;
	model.gradient(initX, grad);

	std::cout << "started energy: " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << energy << ", ||g|| = " << grad.template lpNorm<Eigen::Infinity>() << std::endl;

	if(grad.norm() < params.gradNorm)
	{
		std::cout << "init gradient norm is small" << std::endl;
		return;
	}
	cppoptlib::NewtonDescentSolver<ElasticShellModel> newtonsolver;

	cppoptlib::Criteria<double> s;
	s.iterations = params.iterations;
	s.fDelta = params.fDelta;
	s.gradNorm = params.gradNorm;
	s.xDelta = params.xDelta;
	// if(setup.tensionField) 
	// 	s.fDelta = 1e-10;
	newtonsolver.setStopCriteria(s);
	std::cout << "Stop Criteria: " << std::endl;
	std::cout << "Max iterations: " << s.iterations << ", ||g|| = " << s.gradNorm << ", ||delta_f|| = " << s.fDelta << ", ||delta_x|| = " << s.xDelta << std::endl;

	newtonsolver.minimize(model, initX);
	energy = model.value(initX);
	model.gradient(initX, grad);
	std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << "end energy: " << energy << ", ||grad||: " << grad.norm() << std::endl;
	model.convertVariables2CurState(initX, curState);
	model.convertCurState2Variables(curState, initX);

}