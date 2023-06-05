#include  <igl/boundary_loop.h>
#include <igl/writeOBJ.h>
#include <igl/readOBJ.h>
#include "ShellSolver.h"
#include "Optimization/NewtonDescent.h"
#include "Optimization/NASOQSQP.h"


void ShellSolver::fullSimNewtonStaticSolver(const ElasticSetup& setup, ElasticState& curState, std::string filePrefix, const FullSimOptimizationParams& params)
{
	ElasticShellModel model;
	bool ok = model.initialization(setup, curState, filePrefix, params.isProjH, params.isParallel);
	if (!ok)
    {
        std::cout << "initialization failed." << std::endl;
        return;
    }

	Eigen::VectorXd initX;
	Eigen::MatrixXi F;
	model.convertCurState2Variables(curState, initX);
	double energy = model.value(initX);

	Eigen::VectorXd grad;
	model.gradient(initX, grad);

	if(params.printLog)
	{
		std::cout << "started energy: " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << energy << ", ||g|| = " << grad.template lpNorm<Eigen::Infinity>() << std::endl;
	}


	if(grad.norm() < params.gradNorm)
	{
		if(params.printLog)
			std::cout << "init gradient norm is small" << std::endl;
		return;
	}

    auto elasticFunc = [&](const Eigen::VectorXd& x, Eigen::VectorXd* deriv, Eigen::SparseMatrix<double>* hess, bool isProj)
    {
        model._isUsePosHess = isProj;
        double energy = model.value(x);

        if(deriv)
            model.gradient(x, *deriv);
        if (hess)
            model.hessian(x, *hess);
        return energy;
    };

    OptSolver::newtonSolver(elasticFunc, initX, params.iterations, params.gradNorm, params.xDelta, params.fDelta, params.printLog);
    model.convertVariables2CurState(initX, curState);
}

void ShellSolver::TFWSQPSolver(const TFWSetup& setup, TFWState& curState, std::string filePrefix, const TFWOptimizationParams& params)
{
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
		if(params.printLog)
			std::cout << "negative amplitude, use abs value instead!" << std::endl;
		curState.amplitude = curState.amplitude.cwiseAbs();
	}

	TFWModel model;
	model.initialization(setup, curState, &clampedDOFs, filePrefix, true, params.isParallel);

	Eigen::VectorXd initX;
	model.convertCurState2Variables(curState, initX);
	std::cout << initX.norm() << std::endl;

	auto C = model.buildIntegrabilityConstraints();
	double energy = model.value(initX);

	if(params.printLog)
	{
		std::cout << "initial residual of integrabilty : " << (C * initX).norm() << std::endl;
		std::cout << "started energy: " << energy << std::endl;
	}


	auto elasticFunc = [&](const Eigen::VectorXd& x, Eigen::VectorXd* deriv, Eigen::SparseMatrix<double>* hess, bool isProj)
	{
		model._isUsePosHess = isProj;
		double energy = model.value(x);

		if(deriv)
			model.gradient(x, *deriv);
		if (hess)
			model.hessian(x, *hess);
		return energy;
	};

	auto buildConstraints = [&](const Eigen::VectorXd& x, Eigen::SparseMatrix<double>& Aeq, Eigen::VectorXd& Beq, Eigen::SparseMatrix<double>& Aineq, Eigen::VectorXd& Bineq)
	{
		Aeq = C;
		int nFreeAmps = model.nFreeAmp();
		std::vector<Eigen::Triplet<double>> AineqCoeff;
		AineqCoeff.clear();
		for (int i = 0; i < nFreeAmps; i++)
		{
			AineqCoeff.emplace_back(i, i, -1.0);    // amp >= 0 <=> -amp <= 0
		}
		Aineq.resize(nFreeAmps, x.size());
		Aineq.setFromTriplets(AineqCoeff.begin(), AineqCoeff.end());

		Beq = Eigen::VectorXd::Zero(Aeq.rows());
		Bineq = Eigen::VectorXd::Zero(Aineq.rows());
		return;
	};

	OptSolver::NASOQ_SQPSolver(elasticFunc, initX, buildConstraints, params.iterations, params.gradNorm, params.xDelta, params.fDelta, params.printLog);

	if(params.printLog)
	{
		std::cout << "final residual of integrabilty : " << (C * initX).norm() << std::endl;
	}
	model.convertVariables2CurState(initX, curState);
}
