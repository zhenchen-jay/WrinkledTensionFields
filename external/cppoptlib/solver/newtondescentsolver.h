// CppNumericalSolver
#ifndef NEWTONDESCENTSOLVER_H_
#define NEWTONDESCENTSOLVER_H_

#include <iostream>
#include <Eigen/LU>
#include <Eigen/CholmodSupport>
#include "isolver.h"
#include "../linesearch/morethuente.h"
#include "../linesearch/armijo.h"

#include "../timer.h"
#include "../../../CommonFunctions.h"

namespace cppoptlib {
template<typename ProblemType>
class NewtonDescentSolver : public ISolver<ProblemType, 2> {
  public:
	using Superclass = ISolver<ProblemType, 2>;
	using typename Superclass::Scalar;
	using typename Superclass::TVector;
	using typename Superclass::THessian;
	void minimize(ProblemType& objFunc, TVector& x0) {
		TimeCost totalTimingRecord;
		const int DIM = x0.rows();
		TVector grad = TVector::Zero(DIM);
		THessian hessian;
		TVector neggrad, delta_x;
		Scalar alpha_init = 1.0;
		double maxStepSize = 1.0;

		this->m_current.reset();
		double reg = 1e-6;
		objFunc.gradient(x0, grad);
		std::cout << "initial ||g||_inf = " << grad.template lpNorm<Eigen::Infinity>() << std::endl;
		objFunc.save(this->m_current.iterations, totalTimingRecord, 0, 0, 0, 0, 0, 0, objFunc._isUsePosHess);	// save initial state

		do {
			Timer<std::chrono::high_resolution_clock> timer;
			TimeCost timingRecord;
			std::cout << "iter: " << this->m_current.iterations << std::endl;
			bool isPosHess = objFunc._isUsePosHess;
			if (isPosHess)
				std::cout << "Using positive projected hessian." << std::endl;
			else
				std::cout << "Actual hessian." << std::endl;

			timer.start();
			objFunc.gradient(x0, grad);
			timer.stop();
			timingRecord.gradTime = timer.elapsed<std::chrono::milliseconds>() * 1e-3;

			std::cout << "gradient took: " << timingRecord.gradTime << std::endl;

			timer.start();
			objFunc.hessian(x0, hessian);
			timer.stop();
			timingRecord.hessTime = timer.elapsed<std::chrono::milliseconds>() * 1e-3;
			std::cout << "hessian took: " << timingRecord.hessTime << std::endl;

			timer.start();

			Eigen::SparseMatrix<double> H = hessian;
			Eigen::SparseMatrix<double> I(DIM, DIM);
			I.setIdentity();
			hessian = H + reg * I;
			timer.start();
			Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double> > solver(hessian);
			// Eigen::CholmodSimplicialLLT<Eigen::SparseMatrix<double> > solver(hessian);

			while (solver.info() != Eigen::Success)			// note that we don't PD-project the pressure hessian
			{
				reg = std::max(2 * reg, 1e-16);
				std::cout << "Matrix is not positive definite, current reg = " << reg << std::endl;
				hessian = H + reg * I;
				solver.compute(hessian);
			}
			timer.stop();
			timingRecord.solverTime = timer.elapsed<std::chrono::milliseconds>() * 1e-3;
			std::cout << "LLT took: " << timingRecord.solverTime << std::endl;

			neggrad = -grad;
			delta_x = solver.solve(neggrad);
			alpha_init = std::min(1.0 / grad.norm(), 1.0);
			double descent = grad.dot(delta_x);
			while (descent > 0)
			{
				if(reg < 1e4)
				{
					reg *= 2;
					std::cout << "Not a descent direction. Increase reg to make H more PD!" << std::endl;
					hessian = H + reg * I;
					solver.compute(hessian);
					delta_x = solver.solve(neggrad);
					descent = grad.dot(delta_x);
					alpha_init = std::min(1.0 / grad.norm(), 1.0);
				}
				else
				{
					std::cout << "Not a descent directio and reg is too high > 1e4. Use gradient descent instead!" << std::endl;
					alpha_init = 1.0;
					delta_x = -grad;
				}
			}

			timer.start();
			maxStepSize = objFunc.getMaxStep(x0, delta_x, 1.0);
			timer.stop();
			timingRecord.collisionDectionTime = timer.elapsed<std::chrono::milliseconds>() * 1e-3;
			std::cout << "find max step size took: " << timingRecord.collisionDectionTime << ", maximum time step = " << maxStepSize << ", g^T d = " << descent << std::endl;


			if (maxStepSize > 1e-5)
			{
				std::cout << "current penalty stiffness is : " << objFunc._setup.penaltyK << std::endl;
			}
			else
			{
				delta_x = -grad;
				timer.start();
				maxStepSize = objFunc.getMaxStep(x0, delta_x, 1.0);
				timer.stop();
				std::cout << "last attempt with gradient descent! " << std::endl;
				timingRecord.collisionDectionTime += timer.elapsed<std::chrono::milliseconds>() * 1e-3;
				std::cout << "find max step size took: " << timingRecord.collisionDectionTime << ", maximum time step = " << maxStepSize << std::endl;
			}
			//}


			// find steplength
			Scalar rate = 0;
			// if(!objFunc.isC2())
			//     rate = Armijo<ProblemType, 1>::linesearch(x0, delta_x, objFunc, alpha_init);
			// else
			timer.start();
			rate = MoreThuente<ProblemType, 1>::linesearch(x0, delta_x, objFunc, alpha_init, maxStepSize);

			if(rate == 0 && descent <= 0)
				rate = Armijo<ProblemType, 1>::linesearch(x0, delta_x, objFunc, alpha_init);

			else if (rate == 0)
			{
				std::cout << "Not a descent direction. Use gradient descent instead!" << std::endl;
				alpha_init = 1.0;
				delta_x = -grad;
				// if(!objFunc.isC2())
				//     rate = Armijo<ProblemType, 1>::linesearch(x0, delta_x, objFunc, alpha_init);
				// else
				rate = MoreThuente<ProblemType, 1>::linesearch(x0, delta_x, objFunc, alpha_init, maxStepSize);
				// rate = Armijo<ProblemType, 1>::linesearch(x0, delta_x, objFunc, alpha_init);
			}

			timer.stop();
			timingRecord.lineSearchTime = timer.elapsed<std::chrono::milliseconds>() * 1e-3;
			std::cout << "line search took: " << timingRecord.lineSearchTime << std::endl;

			timer.start();
			reg *= 0.5;
			reg = std::max(reg, 1e-16);

			TVector x1 = x0;
			x0 = x0 + rate * delta_x;

			timer.stop();
			timingRecord.updateTime = timer.elapsed<std::chrono::milliseconds>() * 1e-3;
			std::cout << "update took: " << timingRecord.updateTime << std::endl;


			timer.start();
			double f_old = objFunc.value(x1);
			double f_new = objFunc.value(x0);
			objFunc.gradient(x0, grad);

			std::cout << std::endl << "line search rate = " << rate << ", initial rate = " << alpha_init << ", max step size = " << maxStepSize << ", reg " << reg << std::endl;
			std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << "f_old = " << f_old << ", f_new = " << f_new << std::setprecision(10) << ", delta_f: " << f_old - f_new << ", ||g||:  " << grad.template lpNorm<Eigen::Infinity>() << ", ||dir||: " << delta_x.template lpNorm<Eigen::Infinity>() << std::endl;
			double stretchingEnergy = objFunc.stretchingValue(x0);
			double bendingEnergy = objFunc.bendingValue(x0);
			std::cout << "stretching energy: " << stretchingEnergy << ", bending energy: " << bendingEnergy << std::endl << std::endl;

			++this->m_current.iterations;
			this->m_current.gradNorm = grad.template lpNorm<Eigen::Infinity>();
			this->m_current.xDelta = (rate * delta_x).template lpNorm<Eigen::Infinity>();
			this->m_current.fDelta = f_old - f_new;
			this->m_status = checkConvergence(this->m_stop, this->m_current);

			timer.stop();
			timingRecord.convergenceCheckTime = timer.elapsed<std::chrono::milliseconds>() * 1e-3;
			std::cout << "convergence check took: " << timingRecord.convergenceCheckTime << std::endl;
			
			totalTimingRecord.gradTime += timingRecord.gradTime;
			totalTimingRecord.hessTime += timingRecord.hessTime;
			totalTimingRecord.solverTime += timingRecord.solverTime;
			totalTimingRecord.collisionDectionTime += timingRecord.collisionDectionTime;
			totalTimingRecord.lineSearchTime += timingRecord.lineSearchTime;
			totalTimingRecord.convergenceCheckTime += timingRecord.convergenceCheckTime;

			objFunc.save(this->m_current.iterations, totalTimingRecord, rate, f_old, f_new, grad.template lpNorm<Eigen::Infinity>(), delta_x.template lpNorm<Eigen::Infinity>(), reg, objFunc._isUsePosHess);

			timer.stop();
			timingRecord.savingTime = timer.elapsed<std::chrono::milliseconds>() * 1e-3;
			std::cout << "saving took: " << timingRecord.savingTime << std::endl;

			
			// totalTimingRecord.savingTime += timingRecord.savingTime;

			std::cout << "total time spent: " << totalTimingRecord.totalTime() << std::endl;

			if ((grad.template lpNorm<Eigen::Infinity>() < 1e-7 || rate * delta_x.template lpNorm<Eigen::Infinity>() < 1e-7 || this->m_current.fDelta < 1e-5) && objFunc._isUsePosHess)
			{
				std::cout << "Near the optimal, switch to the actual hessian." << std::endl;
				objFunc._isUsePosHess = false;
			}

		} while (objFunc.callback(this->m_current, x0) && (this->m_status == Status::Continue));
	}
};
}
/* namespace cppoptlib */
#endif /* NEWTONDESCENTSOLVER_H_ */