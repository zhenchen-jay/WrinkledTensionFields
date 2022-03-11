// CppNumericalSolver
#pragma once

#include <iostream>
#include <Eigen/LU>
#include <Eigen/CholmodSupport>
#include "isolver.h"
#include "../linesearch/morethuente.h"
#include "../linesearch/armijo.h"
#include "../../../EigenNASOQ.h"
#include "../timer.h"
#include "../../../CommonFunctions.h"



namespace cppoptlib {
	template<typename ProblemType>
	class NASOQSolver : public ISolver<ProblemType, 2> {
	public:
		using Superclass = ISolver<ProblemType, 2>;
		using typename Superclass::Scalar;
		using typename Superclass::TVector;
		using typename Superclass::THessian;

		void minimize(ProblemType& objFunc, TVector& x0) {
			Timer<std::chrono::high_resolution_clock> totalTimer;
			TimeCost totalTimingRecord;

			totalTimer.start();

			const int DIM = x0.rows();
			TVector grad = TVector::Zero(DIM);
			THessian hessian;
			this->m_current.reset();
			this->m_status = checkConvergence(this->m_stop, this->m_current);

			double reg = 1e-6;
			double PDreg = 1e-6;
			double SPDH_diagPerturb = 1e-9;
			double actualH_diagPerturb = 1e-9;
			double nasoq_eps = std::min(1e-6, this->m_stop.gradNorm);
			std::cout << "nasoq eps: " << nasoq_eps << std::endl;

			int nAmps = objFunc.nFreeAmp();

			Eigen::VectorXd B, Beq, Bieq;
			Eigen::SparseMatrix<double> Q, Aeq, Aieq, C;
			Eigen::VectorXd delta_x;

			int nVars = x0.size();
			int nEq = 0, nIneq = nAmps;
			Aieq.resize(nIneq, nVars);
			Bieq.resize(nIneq);

			// form the selection matrix
			std::vector<Eigen::Triplet<double>> AineqCoeff;
			AineqCoeff.clear();
			for (int i = 0; i < nIneq; i++)
			{
				AineqCoeff.emplace_back(i, i, -1.0);
			}
			Aieq.setFromTriplets(AineqCoeff.begin(), AineqCoeff.end());
			delta_x.resize(nVars);
			delta_x.setZero();

			double totalSavingTime = 0;
			double actualSolverTime = 0;

			C = objFunc.buildIntegrabilityConstraints();
			nEq = C.rows();
			std::cout << "number of equality constraints: " << nEq << std::endl;

			objFunc._isUsePosHess = true;
			bool isUseGradient = false;
			bool notDescentHappened = false;

			Eigen::VectorXd dualy(nEq), dualz(nIneq);
			dualy.setZero();
			dualz.setZero();

			// save initial state
			objFunc.save(0, totalTimingRecord, 0, objFunc.value(x0), objFunc.value(x0), 0, 0, 0, 0, 0, objFunc._isUsePosHess);

			do {
				Timer<std::chrono::high_resolution_clock> timer, localTimer;
				TimeCost timingRecord;

				std::cout << std::endl << "iter: " << this->m_current.iterations << std::endl;

				localTimer.start();
				timer.start();
				objFunc.gradient(x0, grad);
				timer.stop();
				timingRecord.gradTime = timer.elapsed<std::chrono::milliseconds>() * 1e-3;
				totalTimingRecord.gradTime += timingRecord.gradTime;

				std::cout << "gradient took: " << timingRecord.gradTime << std::endl;

				timer.start();
				objFunc.hessian(x0, hessian);
				timer.stop();
				timingRecord.hessTime = timer.elapsed<std::chrono::milliseconds>() * 1e-3;
				totalTimingRecord.hessTime += timingRecord.hessTime;
				std::cout << "hessian took: " << timingRecord.hessTime << std::endl;
				std::cout << "hessian norm: " << hessian.norm() << std::endl;

				timer.start();

				Eigen::SparseMatrix<double> H = hessian;

				Eigen::SparseMatrix<double> I(DIM, DIM);
				I.setIdentity();

				if (isUseGradient)
				{
					std::cout << "use gradient descent, setting hessian to identity!" << std::endl;
					hessian = I;
				}
				else if (!objFunc._isUsePosHess)
				{
					std::cout << "Use actual hessian." << std::endl;
					hessian = H + reg * I;
					Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double> > solver(hessian);
					while (solver.info() != Eigen::Success)
					{
						reg = std::max(2 * reg, 1e-16);
						std::cout << "Matrix is not positive definite, current reg = " << reg << std::endl;
						hessian = H + reg * I;
						solver.compute(hessian);
					}
				}
				else
				{
					std::cout << "Use the projected local hessian." << std::endl;
				}
				Scalar rate = 0;

				std::cout << "minimal amp value: " << x0.segment(0, nAmps).minCoeff() << std::endl;

				Aeq = C;
				nEq = Aeq.rows();
				std::cout << "Integralbility check: " << (C * x0).norm() << std::endl;
				Beq = -C * x0;

				for (int i = 0; i < nAmps; i++)		// x_a + d_a >= 0 <=> S(x + d) >= 0 <=> S*d >= -S*x = -x_a <=> (-S) * d <= x_a
					Bieq(i) = x0(i);

				if (notDescentHappened)
				{
					for (int i = 0; i < nAmps; i++)
						Bieq(i) = x0(i) > 0 ? x0(i) : 0;
					Beq.setZero();
				}


				delta_x = Eigen::VectorXd::Zero(x0.size(), 1);
				Q = hessian;
				B = grad.sparseView();

				double value = -1;
				double QPTime = 0;
				double lineSearchTime = 0;

				std::cout << "QP Solver Start: " << std::endl;
				EigenNASOQSparse qp;
				if (notDescentHappened)
				{
					std::cout << "pure QP. " << std::endl;
					qp.setAccThresh(nasoq_eps * 0.1);
				}
				else
					qp.setAccThresh(nasoq_eps);
				bool isQPSuccess = true;
				Eigen::VectorXd sqpDualy, sqpDualz;

				if (objFunc._isUsePosHess)	// make sure Q is positive definite, not just PD
				{
					Q = Q + PDreg * I;
					SPDH_diagPerturb = std::min(SPDH_diagPerturb, 1e-6);
					isQPSuccess = qp.solve(Q, B, Aeq, Beq, Aieq, Bieq, SPDH_diagPerturb, delta_x, sqpDualy, sqpDualz);
				}
				else
				{
					actualH_diagPerturb = std::min(actualH_diagPerturb, 1e-6);
					isQPSuccess = qp.solve(Q, B, Aeq, Beq, Aieq, Bieq, actualH_diagPerturb, delta_x, sqpDualy, sqpDualz);
				}


				timer.stop();
				timingRecord.solverTime = timer.elapsed<std::chrono::milliseconds>() * 1e-3;
				std::cout << "QP solver took: " << timingRecord.solverTime << std::endl;
				totalTimingRecord.solverTime += timingRecord.solverTime;

				double descent = B.dot(delta_x);
				value = std::abs(descent);
				std::cout << "0.5 * d^T * H * d + g^T * d = " << 0.5 * delta_x.transpose() * Q * delta_x + B.dot(delta_x) << std::endl;
				std::cout << "0.5 * d^T * H * d = " << 0.5 * delta_x.transpose() * Q * delta_x << std::endl;
				std::cout << "g^T * d = " << descent << std::endl;

				std::cout << "(amp + delta_a)_min = " << (delta_x + x0).segment(0, nAmps).minCoeff() << std::endl;
				if (Aeq.rows() > 0)
				{
					std::cout << "||C * delta_x - Beq|| = " << (C * (delta_x + x0)).norm() << std::endl;
				}
				if (!isQPSuccess)
				{
					if (isUseGradient)	// just gradient descent
					{
						std::cout << "QP solver failed!" << std::endl;
						break;
					}
					if (objFunc._isUsePosHess)	// PD projection case
					{
						std::cout << "SPD hessian failed, increase the PD reg, current reg = " << PDreg << ", ||H|| = " << H.norm() << std::endl;
						if (PDreg > 1e4)
						{
							std::cout << "reg is too huge, switch to gradient descent" << std::endl;
							isUseGradient = true;
							localTimer.stop();
							actualSolverTime += localTimer.elapsed<std::chrono::milliseconds>() * 1e-3;
							continue;
						}
						PDreg = std::max(2 * PDreg, 1e-16);
						localTimer.stop();
						actualSolverTime += localTimer.elapsed<std::chrono::milliseconds>() * 1e-3;
						continue;
					}
					else		// actual hessian case
					{
						std::cout << "actual hessian with QP solver failure, increase reg, current reg = " << reg << ", ||H|| = " << H.norm() << std::endl;
						if (reg > 1e4)
						{
							std::cout << "reg is too huge, switch to gradient descent" << std::endl;
							isUseGradient = true;
							localTimer.stop();
							actualSolverTime += localTimer.elapsed<std::chrono::milliseconds>() * 1e-3;
							continue;
						}
						reg = std::max(2 * reg, 1e-16);
						localTimer.stop();
						actualSolverTime += localTimer.elapsed<std::chrono::milliseconds>() * 1e-3;
						continue;
					}

				}
				else if (descent > 0) // QP succeed but, not a descent direction
				{
					if (notDescentHappened && isUseGradient)	// already tried with pure QP and used gradient descent
					{
						std::cout << "solver failed" << std::endl;
						break;
					}
					else if (notDescentHappened)
					{
						std::cout << "pure QP with hessian failed, last attempt: gradient descent" << std::endl;
						isUseGradient = true;
						localTimer.stop();
						actualSolverTime += localTimer.elapsed<std::chrono::milliseconds>() * 1e-3;
						continue;
					}
					else
					{
						std::cout << "Not a descent direction. try with pure QP." << std::endl;
						localTimer.stop();
						actualSolverTime += localTimer.elapsed<std::chrono::milliseconds>() * 1e-3;
						/*
						previous, we solve min_d 1/2 d^T H d + g^Td
										s.t.     d_a + x_a >= 0
												A(d + x) = 0
						By KKT condition, we have,
									g^T d = - d^T H d + \lambda^T A x0 + \mu^T x_a
										= - d^T H d + \mu^T x_a < -d^T H d <= 0 (\mu <= 0)
						But, due to the tolerance in numerical, this might be false near the optimal (A x0 != 0, and x_a has some negative component). Once this happened, we switcch to solve
											min_d 1/2 d^T H d + g^Td
										s.t.     d_a + max(0, x_a) >= 0
												A d = 0
						this will guarantee the descent.
						*/

						notDescentHappened = true;
						continue;
					}
				}
				else
				{
					std::cout << "QP succeeded to find a descent direction. " << std::endl;
					// notDescentHappened = false;
					timer.start();
					// find steplength
					rate = MoreThuente<ProblemType, 1>::linesearch(x0, delta_x, objFunc, 1.0, 1.0);
					timer.stop();

					timingRecord.lineSearchTime = timer.elapsed<std::chrono::milliseconds>() * 1e-3;
					totalTimingRecord.lineSearchTime += timingRecord.lineSearchTime;
					std::cout << "line search took: " << timingRecord.lineSearchTime << std::endl;


					timer.start();
					if (rate <= 1e-10) // due to the tolerance error, the direction may be blurred.
					{
						std::cout << "line search failed with step size =" << rate << std::endl;
						notDescentHappened = true;
						continue;
					}
					else
						isUseGradient = false;
				}

				// Remark: this part can be removed if don't want to record the unconstrained LLT timing
				timer.start();
				if (objFunc._isUsePosHess)
					Q = Q + 1e-8 * I;
				Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double> > solver(Q);

				Eigen::VectorXd negGrad = -grad;
				Eigen::VectorXd freeDir = solver.solve(negGrad);
				timer.stop();
				std::cout << "LLT unconstrained solver took: " << timer.elapsed<std::chrono::milliseconds>() * 1e-3 << std::endl;
				timingRecord.unconstrainedLLTTime = timer.elapsed<std::chrono::milliseconds>() * 1e-3;
				totalTimingRecord.unconstrainedLLTTime += timingRecord.unconstrainedLLTTime;



				TVector x1 = x0;
				// update primal variables
				x0 = x0 + rate * delta_x;
				// update dual variables
				dualy = dualy + rate * (sqpDualy - dualy);
				dualz = dualz + rate * (sqpDualz - dualz);

				timer.stop();
				timingRecord.updateTime = timer.elapsed<std::chrono::milliseconds>() * 1e-3;
				std::cout << "update took: " << timingRecord.updateTime << std::endl;

				timer.start();
				double f_old = objFunc.value(x1);
				double f_new = objFunc.value(x0);

				double f_delta = f_old - f_new;
				double x_delta = (x1 - x0).template lpNorm<Eigen::Infinity>();

				double amp_delta = (x1 - x0).segment(0, objFunc.nFreeAmp()).template lpNorm<Eigen::Infinity>();
				double dphi_delta = (x1 - x0).segment(objFunc.nFreeAmp(), x1.size() - objFunc.nFreeAmp()).template lpNorm<Eigen::Infinity>();

				double amp_inf = x0.segment(0, objFunc.nFreeAmp()).template lpNorm<Eigen::Infinity>();
				double dphi_inf = x0.segment(objFunc.nFreeAmp(), x1.size() - objFunc.nFreeAmp()).template lpNorm<Eigen::Infinity>();

				TVector curGrad;
				double subGradInfNorm = std::numeric_limits<double>::infinity();
				int nSize = x0.size();
				objFunc.gradient(x0, curGrad);

				std::cout << std::endl << "line search rate = " << rate;
				if (objFunc._isUsePosHess)
					std::cout << ", PD reg " << PDreg << std::endl;
				else
					std::cout << ", reg " << reg << std::endl;
				std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << "f_old = " << f_old << ", f_new = " << f_new << std::setprecision(10) << std::endl;
				std::cout << "|<grad, dir>| = " << value << ", f_delta = " << f_delta << ", x_delta = " << x_delta << ", amp_delta = " << amp_delta << ", dphi_delta = " << dphi_delta << std::endl;

				double stretchingEnergy = objFunc.stretchingValue(x0);
				double bendingEnergy = objFunc.bendingValue(x0);
				std::cout << "stretching energy: " << stretchingEnergy << ", bending energy: " << bendingEnergy << std::endl;

				std::cout << "||amp||_inf : " << amp_inf << ", ||dphi||_inf: " << dphi_inf << std::endl << std::endl;

				std::cout << "lagrangian informations: " << std::endl;
				subGradInfNorm = (curGrad + Aeq.transpose() * dualy + Aieq.transpose() * dualz).template lpNorm<Eigen::Infinity>();
				std::cout << "stationarity check: " << subGradInfNorm << std::endl;
				Eigen::VectorXd ineqCheck = Aieq * x0;
				for (int i = 0; i < nIneq; i++)
					ineqCheck(i) = ineqCheck(i) > 0 ? ineqCheck(i) : 0;
				Eigen::VectorXd dualCheck = dualz;
				for (int i = 0; i < dualz.size(); i++)
					dualCheck(i) = dualz(i) < 0 ? dualz(i) : 0;
				std::cout << "primal feasibility: " << std::endl << "Equality part: " << (Aeq * x0).template lpNorm<Eigen::Infinity>() << ", Inequality part: " << ineqCheck.maxCoeff() << std::endl;
				std::cout << "dual feasibilty: " << dualCheck.template lpNorm<Eigen::Infinity>() << std::endl;
				std::cout << "complemtarity: " << (dualz.cwiseProduct(Aieq * x0)).template lpNorm<Eigen::Infinity>() << std::endl;



				++this->m_current.iterations;
				this->m_current.gradNorm = subGradInfNorm;
				this->m_current.xDelta = x_delta;
				this->m_current.fDelta = f_delta;
				this->m_status = checkConvergence(this->m_stop, this->m_current);

				timer.stop();
				timingRecord.convergenceCheckTime = timer.elapsed<std::chrono::milliseconds>() * 1e-3;
				totalTimingRecord.convergenceCheckTime += timingRecord.convergenceCheckTime;
				std::cout << "convergence check took: " << timingRecord.convergenceCheckTime << std::endl;

				localTimer.stop();
				actualSolverTime += localTimer.elapsed<std::chrono::milliseconds>() * 1e-3;

				timer.start();
				if (objFunc._isUsePosHess)
					objFunc.save(this->m_current.iterations, totalTimingRecord, rate, f_old, f_new, f_delta, amp_delta, dphi_delta, PDreg, &subGradInfNorm, objFunc._isUsePosHess, true); // save all the intermediate results
				else
					objFunc.save(this->m_current.iterations, totalTimingRecord, rate, f_old, f_new, f_delta, amp_delta, dphi_delta, reg, &subGradInfNorm, objFunc._isUsePosHess, true); // save all the intermediate results
				timer.stop();
				double savingTime = timer.elapsed<std::chrono::milliseconds>() * 1e-3;
				std::cout << "saving took: " << savingTime << std::endl;
				totalSavingTime += savingTime;

				if (!isUseGradient)
				{
					if (objFunc._isUsePosHess)
					{
						PDreg *= 0.5;
						PDreg = std::max(PDreg, 1e-16);
					}
					else
					{
						reg *= 0.5;
						reg = std::max(reg, 1e-16);
					}
				}


				if (f_delta < 1e-8) // near the local minima and the actual hessian is basically PD
				{
					objFunc._isUsePosHess = false;
				}
				else
				{
					objFunc._isUsePosHess = true;
				}
				std::cout << "total time spent: " << totalTimer.elapsed<std::chrono::milliseconds>() * 1e-3 << ", actual solver time: " << actualSolverTime << ", saving time: " << totalSavingTime << std::endl;

			} while (objFunc.callback(this->m_current, x0) && (this->m_status == Status::Continue));

			for (int i = 0; i < objFunc.nFreeAmp(); i++)
			{
				if (x0(i) < 0)
					x0(i) = 0;
			}
		}
	};
}
