#include <Eigen/CholmodSupport>
#include <fstream>
#include <iomanip>
#include "LineSearch.h"
#include "NASOQSQP.h"
#include "../EigenNASOQ.h"
#include "../Timer.h"

void OptSolver::NASOQ_SQPSolver(
		std::function<double(const Eigen::VectorXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool)> objFunc, // objective function
		Eigen::VectorXd& x0, // start point (also the optimal point)
		std::function<void(const Eigen::VectorXd&, Eigen::SparseMatrix<double>&, Eigen::VectorXd&, Eigen::SparseMatrix<double>&, Eigen::VectorXd&)> constraintFunc, // constraint function: constraintFunc(x, Aeq, Beq, Aineq, Bineq);
		int numIter,     // maximum iteration
		double gradTol, // gradient norm tolerance
		double xTol,        // variable update tolerance
		double fTol,        // function value update tolerance
		bool displayInfo,   // whether display the optimization information every step
		std::function<void(const Eigen::VectorXd&, double&, double&)> getNormFunc,    // optional function for get the norm
		std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&)> findMaxStep,    // optional function for find the maximum step, like collision detection
		std::string *savingFolder,        // whether save the intermediate results. If provided, save all intermediate results to the that folder
		std::function<void(const Eigen::VectorXd&, std::string*)> saveProcess     // optional function for save the problem
)
{
	Timer totalTimer;
	double totalAssemblingTime = 0;
	double totalSolvingTime = 0;
	double totalLineSearchTime = 0;
	double totalConvergenceCheckTime = 0;
	totalTimer.start();
	std::ofstream optInfo;

	auto coutLog = [&](const std::string& str)
	{
		if(displayInfo)
		{
			std::cout << str << std::endl;
		}
		if(saveProcess)
		{
			optInfo << str << std::endl;
		}
	};

	auto double2Str = [&](const double& val, unsigned int numDigits = 6)
	{
		std::ostringstream streamObj;
		streamObj << std::fixed;
		streamObj << std::setprecision(numDigits);
		streamObj << val;
		return streamObj.str();
	};

	if (savingFolder)
	{
		optInfo = std::ofstream((*savingFolder) + "optInfo.txt");
	}
	coutLog("NASOQ SQP solver with termination criterion: ");
	coutLog("gradient tol: " + std::to_string(gradTol) + ", function update tol: " + std::to_string(fTol) + ", variable update tol: " + std::to_string(xTol) + ", maximum iteration: " + std::to_string(numIter));

	Eigen::VectorXd grad;
	Eigen::SparseMatrix<double> Aeq, Aineq, hessian;
	Eigen::VectorXd Beq0, Bineq0;
	constraintFunc(x0, Aeq, Beq0, Aineq, Bineq0);

	const int DIM = x0.size();
	Eigen::VectorXd delta_x = Eigen::VectorXd::Zero(DIM);

	bool isProj = true;
	bool isUseGradient = false;
	bool pureQP = false;

	double reg = 1e-6;
	double PDreg = 1e-6;
	double SPDH_diagPerturb = 1e-9;
	double actualH_diagPerturb = 1e-9;
	double nasoq_eps = std::min(1e-6, gradTol);
	int nEq = Aeq.rows();
	int nIneq = Aineq.rows();

	coutLog("nasoq eps: " + std::to_string(nasoq_eps));
	coutLog("number of equality constraints: " + std::to_string(Aeq.rows()));
	coutLog("number of inequality constraints: " + std::to_string(Aineq.rows()));

	Eigen::VectorXd dualy(nEq), dualz(nIneq);
	dualy.setZero();
	dualz.setZero();

	// given we don't asked for the dual variable as input, we just optimize as it
	int iter = 0;
	Timer localTimer;
	for(iter = 0; iter < numIter; iter++)
	{
		localTimer.start();
		coutLog("\niter: " + std::to_string(iter));
		double f = objFunc(x0, &grad, &hessian, isProj);
		localTimer.stop();
		double localAssTime = localTimer.elapsedSeconds();
		totalAssemblingTime += localAssTime;

		localTimer.start();
		Eigen::SparseMatrix<double> H = hessian;

		Eigen::SparseMatrix<double> I(DIM, DIM);
		I.setIdentity();

		if (isUseGradient)
		{
			if(displayInfo)
				std::cout << "use gradient descent, setting hessian to identity!" << std::endl;
			hessian = I;
		}
		else if (!isProj)
		{
			if(displayInfo)
			{
				std::cout << "Use actual hessian." << std::endl;
			}
			hessian = H + reg * I;
			Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double> > solver(hessian);
			while (solver.info() != Eigen::Success)
			{
				reg = std::max(2 * reg, 1e-16);
				if(displayInfo)
				{
					std::cout << "Matrix is not positive definite, current reg = " << reg << std::endl;
				}
				hessian = H + reg * I;
				solver.compute(hessian);
			}
		}
		else
		{
			if(displayInfo)
				std::cout << "Use the projected local hessian." << std::endl;
		}


		double rate = 0;
		Eigen::VectorXd Beq = Beq0 - Aeq * x0;      // we want Aeq * (x0 + delta_x) = Beq0, since near the optimal, one should expect a unit line search rate
		Eigen::VectorXd Bineq = Beq0 - Aineq * x0;  // again we want Aineq * (x0 + delta_x) <= Beq

		if(pureQP)
		{
			Beq.setZero();

			for (int i = 0; i < Bineq.rows(); i++)
				Bineq(i) = std::max(0., Bineq(i));
		}

		Eigen::SparseMatrix<double> Q;
		Eigen::SparseVector<double> B;

		Q = hessian;
		B = grad.sparseView();
		delta_x.setZero();
		EigenNASOQSparse qp;

		if (pureQP)
		{
			if(displayInfo)
				std::cout << "pure QP. " << std::endl;
			qp.setAccThresh(nasoq_eps * 0.1);
		}
		else
			qp.setAccThresh(nasoq_eps);

		Eigen::VectorXd sqpDualy, sqpDualz;
		bool isQPSuccess = true;

		if(displayInfo)
		{
			std::cout << "QP start" << std::endl;
		}
		if (isProj)	// make sure Q is positive definite, not just positive semi-definite
		{
			Q = Q + PDreg * I;
			SPDH_diagPerturb = std::min(SPDH_diagPerturb, 1e-6);
			isQPSuccess = qp.solve(Q, B, Aeq, Beq, Aineq, Bineq, SPDH_diagPerturb, delta_x, sqpDualy, sqpDualz);
		}
		else
		{
			actualH_diagPerturb = std::min(actualH_diagPerturb, 1e-6);
			isQPSuccess = qp.solve(Q, B, Aeq, Beq, Aineq, Bineq, actualH_diagPerturb, delta_x, sqpDualy, sqpDualz);
		}
		localTimer.stop();
		totalSolvingTime += localTimer.elapsedSeconds();

		if(displayInfo)
		{
			std::cout << "QP solver took: " << totalSolvingTime << std::endl;
		}

		localTimer.start();
		double descent = B.dot(delta_x);
		if(!isQPSuccess)
		{
			if (isUseGradient)	// just gradient descent
			{
				if(displayInfo)
					std::cout << "QP solver failed!" << std::endl;
				break;
			}
			if (isProj)	// PD projection case
			{
				if(displayInfo)
					std::cout << "SPD hessian failed, increase the PD reg, current reg = " << PDreg << ", ||H|| = " << H.norm() << std::endl;
				if (PDreg > 1e4)
				{
					std::cout << "reg is too huge, switch to gradient descent" << std::endl;
					isUseGradient = true;
					localTimer.stop();
					totalSolvingTime += localTimer.elapsedSeconds();
					continue;
				}
				PDreg = std::max(2 * PDreg, 1e-16);
				localTimer.stop();
				totalSolvingTime += localTimer.elapsedSeconds();
				continue;
			}
			else 		// actual hessian case
			{
				if(displayInfo)
					std::cout << "actual hessian with QP solver failure, increase reg, current reg = " << reg << ", ||H|| = " << H.norm() << std::endl;
				if (reg > 1e4)
				{
					std::cout << "reg is too huge, switch to gradient descent" << std::endl;
					isUseGradient = true;
					localTimer.stop();
					totalSolvingTime += localTimer.elapsedSeconds();
					continue;
				}
				reg = std::max(2 * reg, 1e-16);
				localTimer.stop();
				totalSolvingTime += localTimer.elapsedSeconds();
				continue;
			}
		}
		else if (descent > 0) // QP succeed but, not a descent direction
		{
			if (isUseGradient && pureQP)	// already tried with pure QP and used gradient descent
			{
				if(displayInfo)
					std::cout << "solver failed" << std::endl;
				break;
			}
			else if(pureQP)
			{
				if(displayInfo)
					std::cout << "pure QP with hessian failed, last attempt: gradient descent" << std::endl;
				isUseGradient = true;
				localTimer.stop();
				totalSolvingTime += localTimer.elapsedSeconds();
				continue;
			}
			else
			{
				std::cout << "Not a descent direction. try with pure QP." << std::endl;
				localTimer.stop();
				totalSolvingTime += localTimer.elapsedSeconds();
				/*
				previous, we solve min_d 1/2 d^T Q d + g^Td
								s.t.    Ai * d <= (Bi - Ai * x0)
										Ae * d == (Be - Ae * x0)
				Notice that d = 0 is a feasible point theoretically (Ai * x0 <= Bi and Ae * x0 = Be), therefore the optimal d1 should have:
				1/2 d1^T Q d1 + g^Td1 <= 1/2 0^T Q 0 + g^T 0 = 0, which implies g^Td1 <= 1/2 d1^T Q d1 < 0 (Q is PD).
				But, due to numerical error, this might be false near the optimal (A x0 != 0, and Bi - Ai * x0 has some slightly negative component). Once this happened, we switch to solve
									min_d 1/2 d^T Q d + g^Td
								s.t.    Ai * d <= max(0, Bi - Ai * x0)
										Ae * d = 0
				this will guarantee the descent. (make sure d = 0 is the feasible point)
				*/

				pureQP = true;
				continue;
			}
		}
		else
		{
			if(displayInfo)
				std::cout << "QP succeeded to find a descent direction. " << std::endl;
			localTimer.stop();
			totalSolvingTime += localTimer.elapsedSeconds();
			localTimer.start();
			double maxStepSize = findMaxStep ? findMaxStep(x0, delta_x) : 1.0;
			rate = LineSearch::backtrackingArmijo(x0, grad, delta_x, objFunc, maxStepSize);
			localTimer.stop();

			totalLineSearchTime += localTimer.elapsedSeconds();

			localTimer.start();
			if (rate <= 1e-10) // due to the tolerance error, the direction may be blurred.
			{
				if(displayInfo)
				{
					std::cout << "line search failed with step size =" << rate << std::endl;
				}
				pureQP = true;
				continue;
			}
		}

		Eigen::VectorXd x1 = x0, gradNew;
		// update primal variables
		x0 = x0 + rate * delta_x;
		// update dual variables
		dualy = dualy + rate * (sqpDualy - dualy);
		dualz = dualz + rate * (sqpDualz - dualz);

		localTimer.start();
		double fNew = objFunc(x0, &gradNew, nullptr, false);

		double deltaF = f - fNew;
		double deltaX = rate * delta_x.norm();

		double stationarity = (gradNew + Aeq.transpose() * dualy + Aineq.transpose() * dualz).template lpNorm<Eigen::Infinity>();
		Eigen::VectorXd ineqCheck = Aineq * x0 - Bineq0;
		for (int i = 0; i < nIneq; i++)
			ineqCheck(i) = ineqCheck(i) > 0 ? ineqCheck(i) : 0;
		Eigen::VectorXd dualCheck = dualz;
		for (int i = 0; i < dualz.size(); i++)
			dualCheck(i) = dualz(i) < 0 ? dualz(i) : 0;

		coutLog("\nline search rate = " + std::to_string(rate));
		coutLog("Hessian info: (proj hessian, actual hessian, gradient descent): " + std::to_string(isProj) + ", " + std::to_string(!isProj) + ", " + std::to_string(isUseGradient));
		coutLog("PD reg: " + std::to_string(PDreg) + ", reg: " + std::to_string(reg));
		coutLog("f_old = " + double2Str(f, std::numeric_limits<double>::digits10 + 1) + ", f_new = " + double2Str(fNew, std::numeric_limits<double>::digits10 + 1) + ", delta_f = " + double2Str(deltaF, std::numeric_limits<double>::digits10 + 1) + ", delta_x = " + double2Str(deltaX, std::numeric_limits<double>::digits10 + 1));
		coutLog("lagrangian informations: ");
		coutLog("stationarity check: " + double2Str(stationarity, std::numeric_limits<double>::digits10 + 1));
		coutLog("primal feasibility: \nEquality part: " + double2Str((Aeq * x0).template lpNorm<Eigen::Infinity>(), std::numeric_limits<double>::digits10 + 1) + ", Inequality part: " + double2Str(ineqCheck.maxCoeff(), std::numeric_limits<double>::digits10 + 1));
		coutLog("dual feasibilty: " + double2Str(dualCheck.template lpNorm<Eigen::Infinity>(), std::numeric_limits<double>::digits10 + 1));
		coutLog("complemtarity: " + double2Str((dualz.cwiseProduct(Aineq * x0 - Bineq0)).template lpNorm<Eigen::Infinity>(), std::numeric_limits<double>::digits10 + 1));


		// This termination is not actually correct, but practically good
		if(stationarity < gradTol)
		{
			coutLog("Terminate with small stationarity: " + double2Str(stationarity));
			break;
		}
		if(deltaX < xTol)
		{
			coutLog("Terminate with small variable update: " + double2Str(deltaX));
			break;
		}
		if(deltaF < fTol)
		{
			coutLog("Terminate with small function value update: " + double2Str(deltaF));
			break;
		}
		localTimer.stop();
		totalConvergenceCheckTime += localTimer.elapsedSeconds();

		if (!isUseGradient)
		{
			if (isProj)
			{
				PDreg *= 0.5;
				PDreg = std::max(PDreg, 1e-16);
			}
			else
			{
				reg *= 0.5;
				reg = std::max(reg, 1e-16);
			}
		}

		// reset
		pureQP = false;
		isUseGradient = false;

		if (deltaF < 1e-8) // near the local minima and the actual hessian is basically PD
		{
			isProj = false;
		}
		else
		{
			isProj = true;
		}
		coutLog("timing info (in total seconds): \nassembling took: " + std::to_string(totalAssemblingTime) + ", QP solver took: " + std::to_string(totalSolvingTime) + ", line search took: " + std::to_string(totalLineSearchTime) + ", convergence check took: " + std::to_string(totalConvergenceCheckTime));
	}

}