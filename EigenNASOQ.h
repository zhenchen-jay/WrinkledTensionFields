#pragma once

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Core>



class EigenNASOQSparse
{
public:
	EigenNASOQSparse() : _accThresh(1e-6) {}
	bool solve(const Eigen::SparseMatrix<double>& Q, const Eigen::VectorXd& C,
		const Eigen::SparseMatrix<double>& Aeq, const Eigen::VectorXd& Beq,
		const Eigen::SparseMatrix<double>& Aineq, const Eigen::VectorXd& Bineq,
		const Eigen::VectorXd& XL, const Eigen::VectorXd& XU,
		Eigen::VectorXd& X, double& diag_perturbation);
	bool solve(const Eigen::SparseMatrix<double>& Q, const Eigen::VectorXd& C,
		const Eigen::SparseMatrix<double>& Aeq, const Eigen::VectorXd& Beq,
		const Eigen::SparseMatrix<double>& Aineq, const Eigen::VectorXd& Bineq,
		double& diag_perturbation,
		Eigen::VectorXd& X, Eigen::VectorXd& dualy, Eigen::VectorXd& dualz);
	void setAccThresh(double tol) 
	{
		_accThresh = std::max(tol, 1e-16);
	}
	void saveState(const Eigen::SparseMatrix<double>& Q, const Eigen::VectorXd& C,
		const Eigen::SparseMatrix<double>& Aeq, const Eigen::VectorXd& Beq,
		const Eigen::SparseMatrix<double>& Aineq, const Eigen::VectorXd& Bineq,
		const Eigen::VectorXd& XL, const Eigen::VectorXd& XU);
	/*
	solve for the QP: 
		0.5 * x^T * Q * x + c^T x
	s.t.
		Aeq x = Beq
		Aineq x <= Bineq
		x >= xl
		x <= xu
	*/

	// TODO: fix the compile issue
	//  Solving Ax=b
	// Inputs:
	//   H  n by n sparse Hessian matrix **lower triangle only** (see
	//     .triangularView<Eigen::Lower>() )
	//   q  n by 1 vector
	// Outputs:
	//   x  n by 1 solution vector
	// Returns nasoq exit flag
	//
	int linear_solve(
		// Pass inputs by copy so we get non-const and casted data
		Eigen::SparseMatrix<double, Eigen::ColMajor, int> A,
		Eigen::Matrix<double, Eigen::Dynamic, 1> b,
		Eigen::Matrix<double, Eigen::Dynamic, 1>& x);
	void printSolverReturn(int convergedFlag);
private:
	void getNumIter(double acc_thresh, int& inner_iter_ref, int& outer_iter_ref);

	void writeSymMatrix(const Eigen::SparseMatrix<double> H, std::string filename);
	void writeMatrix(const Eigen::SparseMatrix<double> A, std::string filename);
	void writeVector(const Eigen::VectorXd v, std::string filename);

	double _accThresh;

};
