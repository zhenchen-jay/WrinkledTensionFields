#ifndef GUROBIMIPWRAPPER_H
#define GUROBIMIPWRAPPER_H

#include <Eigen/Sparse>
#include <Eigen/Core>

void GurobiMIPWrapper(const Eigen::SparseMatrix<double> &constraints,
    const Eigen::SparseMatrix<double> &A,
    Eigen::VectorXd &result,
    const Eigen::VectorXd &rhs,
    const Eigen::VectorXi &toRound,
    double reg);

#endif
