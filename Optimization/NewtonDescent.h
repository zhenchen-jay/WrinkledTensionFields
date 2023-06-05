#pragma once

#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "TestFunctions.h"


namespace OptSolver
{
	void newtonSolver(std::function<double(const Eigen::VectorXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool)> objFunc, Eigen::VectorXd& x0, int numIter = 1000, double gradTol = 1e-14, double xTol = 0, double fTol = 0, bool displayInfo = false, std::function<void(const Eigen::VectorXd&, double&, double&)> getNormFunc = nullptr, std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&)> findMaxStep = nullptr, std::string *savingFolder = nullptr, std::function<void(const Eigen::VectorXd&, std::string*)> saveProcess = nullptr);

}


