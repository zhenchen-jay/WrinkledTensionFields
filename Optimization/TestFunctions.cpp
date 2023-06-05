#include "TestFunctions.h"

void OptSolver::testFuncGradHessian(std::function<double(const Eigen::VectorXd&, Eigen::VectorXd*, Eigen::SparseMatrix<double>*, bool)> objFunc, const Eigen::VectorXd& x0)
{
	Eigen::VectorXd dir = x0;
	dir(0) = 0;
	dir.setRandom();

	Eigen::VectorXd grad;
	Eigen::SparseMatrix<double> H;

	double f = objFunc(x0, &grad, &H, false);
	std::cout << "f: " << f << std::endl;
    if(f == 0)
        return;

	for (int i = 3; i < 10; i++)
	{
		double eps = std::pow(0.1, i);
		Eigen::VectorXd x = x0 + eps * dir;
		Eigen::VectorXd grad1;
		double f1 = objFunc(x, &grad1, NULL, false);

		std::cout << "\neps: " << eps << std::endl;
		std::cout << "energy-gradient: " << (f1 - f) / eps - grad.dot(dir) << std::endl;
		std::cout << "gradient-hessian: " << ((grad1 - grad) / eps - H * dir).norm() << std::endl;

//		std::cout << "gradient-difference: \n" << (grad1 - grad) / eps << std::endl;
//		std::cout << "direction-hessian: \n" << H * dir << std::endl;
	}
}

