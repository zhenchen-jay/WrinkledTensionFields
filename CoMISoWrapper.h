#ifndef COMISOWRAPPER_H
#define COMISOWRAPPER_H

#include <Eigen/Sparse>
#include <Eigen/Core>

void ComisoWrapper(const Eigen::SparseMatrix<double> &constraints,
    const Eigen::SparseMatrix<double> &A,
    Eigen::VectorXd &result,
    const Eigen::VectorXd &rhs,
    const Eigen::VectorXi &toRound,
    double reg);

#endif
