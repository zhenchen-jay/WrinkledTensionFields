#pragma once

#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "TestFunctions.h"


namespace OptSolver
{
	/* use the Seqeuatial quadratic programming to solve the problem:
	 *      min_x f(x)
	 * s.t.  Aeq x = Beq
	 *       Aineq x <= Bineq
	*/
	void NASOQ_SQPSolver(
			std::function<double(const Eigen::VectorXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool)> objFunc, // objective function
			Eigen::VectorXd& x0, // start point (also the optimal point)
			std::function<void(const Eigen::VectorXd&, Eigen::SparseMatrix<double>&, Eigen::VectorXd&, Eigen::SparseMatrix<double>&, Eigen::VectorXd&)> constraintFunc, // constraint function: constraintFunc(x, Aeq, Beq, Aineq, Bineq);
			int numIter = 1000,     // maximum iteration
			double gradTol = 1e-14, // gradient norm tolerance
			double xTol = 0,        // variable update tolerance
			double fTol = 0,        // function value update tolerance
			bool displayInfo = false,   // whether display the optimization information every step
			std::function<void(const Eigen::VectorXd&, double&, double&)> getNormFunc = nullptr,    // optional function for get the norm
			std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&)> findMaxStep = nullptr,    // optional function for find the maximum step, like collision detection
			std::string *savingFolder = nullptr,        // whether save the intermediate results. If provided, save all intermediate results to the that folder
			std::function<void(const Eigen::VectorXd&, std::string*)> saveProcess = nullptr     // optional function for save the problem
					);

}


