#pragma once

#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Core>


namespace OptSolver
{
	void testFuncGradHessian(std::function<double(const Eigen::VectorXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool)> objFunc, const Eigen::VectorXd& x0);
}