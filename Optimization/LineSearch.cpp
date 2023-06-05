
#include "LineSearch.h"
#include <iostream>


double LineSearch::backtrackingArmijo(const Eigen::VectorXd& x, const Eigen::VectorXd& grad, const Eigen::VectorXd& dir, std::function<double(const Eigen::VectorXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool)> objFunc, const double alphaInit)
{
    const double c = 0.2;
    const double rho = 0.5;
    double alpha = alphaInit;

    Eigen::VectorXd xNew = x + alpha * dir;
    double fNew = objFunc(xNew, NULL, NULL, false);
    double f = objFunc(x, NULL, NULL, false);
    const double cache = c * grad.dot(dir);

    while (fNew > f + alpha * cache) {
        alpha *= rho;
        xNew = x + alpha * dir;
        fNew = objFunc(xNew, NULL, NULL, false);
    }

    return alpha;
}
