#include "CoMISoWrapper.h"
#include <iostream>
#include <CoMISo/Solver/ConstrainedSolver.hh>

void ComisoWrapper(const Eigen::SparseMatrix<double> &constraints,
    const Eigen::SparseMatrix<double> &A, 
    Eigen::VectorXd &result,
    const Eigen::VectorXd &rhs,
    const Eigen::VectorXi &toRound,
    double reg)
{
    int n = A.rows();
    assert(n == A.cols());
    assert(n + 1 == constraints.cols());
    std::cout << n << " " << rhs.size() << std::endl;
    assert(n == rhs.size());
    COMISO::ConstrainedSolver solver;

    gmm::col_matrix< gmm::wsvector< double > > Agmm(n,n);
    int nconstraints = constraints.rows();
    gmm::row_matrix< gmm::wsvector< double > > Cgmm(nconstraints, n+1); // constraints

    for (int k=0; k < A.outerSize(); ++k){
        for (Eigen::SparseMatrix<double>::InnerIterator it(A,k); it; ++it){
            int row = it.row();
            int col = it.col();
            Agmm(row, col) += it.value();
        }
    }
    for (int k=0; k < constraints.outerSize(); ++k){
        for (Eigen::SparseMatrix<double>::InnerIterator it(constraints,k); it; ++it){
            int row = it.row();
            int col = it.col();
            Cgmm(row, col) += it.value();
        }
    }

    std::vector<double> rhsv(n,0); 
    for (int i = 0; i < n; i++)
        rhsv[i] = rhs[i];

    std::vector<double> X(n, 0);
    std::vector<int> toRoundv(toRound.size());
    for (int i = 0; i < toRound.size(); i++)
        toRoundv[i] = toRound[i];
    for (int i = 0; i < toRoundv.size(); i++)
    {
        assert(toRoundv[i] >= 0 && toRoundv[i] < n);
    }
    solver.solve(Cgmm, Agmm, X, rhsv, toRoundv, reg, false, false);
    result.resize(n);
    for (int i = 0; i < n; i++)
        result[i] = X[i];
}

