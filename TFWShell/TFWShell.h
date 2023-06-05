#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <memory>

#include "../MeshLib/MeshConnectivity.h"
#include "../MeshLib/GeometryDerivatives.h"
#include "../MeshLib/MeshGeometry.h"

#include "../SecondFundamentalForm/SecondFundamentalFormDiscretization.h"
#include "../SecondFundamentalForm/MidedgeAverageFormulation.h"
#include "../SecondFundamentalForm/MidedgeAngleSinFormulation.h"
#include "../SecondFundamentalForm/MidedgeAngleTanFormulation.h"


#include "../CommonFunctions.h"

#include "TFWSetup.h"
#include "TFWState.h"

namespace TFW
{
	class TFWShell
	{
	public:
		TFWShell(const TFWSetup& setup,
			const TFWState& state,
			bool isPosHess = false,
			bool isParallel = false)
		{
			_setup = setup;
			_state = state;
			_isPosHess = isPosHess;
			_isParallel = isParallel;
		}

		double stretchingEnergy(Eigen::VectorXd* deriv, Eigen::SparseMatrix<double>* hessian);
		double bendingEnergy(Eigen::VectorXd* deriv, Eigen::SparseMatrix<double>* hessian);

		double elasticReducedEnergy(Eigen::VectorXd* deriv, Eigen::SparseMatrix<double>* hessian);

		double computeAmplitudesFromQuad(int faceId, int quadId, Eigen::Vector2d* da, Eigen::Vector3d* gradA, Eigen::Matrix<double, 2, 3>* gradDA, Eigen::Matrix<double, 3, 3>* hessianA, std::vector<Eigen::Matrix<double, 3, 3>>* hessianDA);
		Eigen::Vector2d computeDphi(int faceId, Eigen::Matrix<double, 2, 3>* gradDphi);

		Eigen::Matrix2d computeDaDphiTensor(int faceId, int quadId, std::vector<Eigen::Matrix2d>* deriv, std::vector<Eigen::Matrix<double, 6, 6>>* hessian);
		Eigen::Matrix2d computeDphiDaTensor(int faceId, int quadId, std::vector<Eigen::Matrix2d>* deriv, std::vector<Eigen::Matrix<double, 6, 6>>* hessian);
		Eigen::Matrix2d computeDaDaTensor(int faceId, int quadId, std::vector<Eigen::Matrix2d>* deriv, std::vector<Eigen::Matrix<double, 3, 3>>* hessian);
		Eigen::Matrix2d computeDphiDphiTensor(int faceId, int quadId, std::vector<Eigen::Matrix2d>* deriv, std::vector<Eigen::Matrix<double, 3, 3>>* hessian);


		std::vector<Eigen::Matrix2d> computeStretchingDensityFromQuad(int faceId, int quadId, std::vector<Eigen::MatrixXd >* deriv, std::vector<std::vector<Eigen::MatrixXd > >* hessian);
		std::vector<Eigen::Matrix2d> computeBendingDensityFromQuad(int faceId, int quadId, std::vector<Eigen::MatrixXd >* deriv, std::vector<std::vector<Eigen::MatrixXd > >* hessian);

		double stretchingEnergyPerface(int faceId, Eigen::VectorXd* deriv, Eigen::MatrixXd* hessian);
		double bendingEnergyPerface(int faceId, Eigen::VectorXd* deriv, Eigen::MatrixXd* hessian);


		Eigen::VectorXd computeStretchingEnergyTermByTerm();
		Eigen::VectorXd computeBendingEnergyTermByTerm();

		// test functions
		void testAmplitudesFromQuad(int faceId, int quadId);
		void testDphi(int faceId);
		void testStretchingDensityFromQuad(int faceId, int quadId);
		void testBendingDensityFromQuad(int faceId, int quadId);
		void testStretchingEnergyPerface(int faceId);
		void testBendingEnergyPerface(int faceId);
		void testStretchingEnergy();
		void testBendingEnergy();
		void testElasticReducedEnergy();


	public:
		TFWSetup _setup;
		TFWState _state;
		MeshGeometry _meshGeo;
		bool _isPosHess;
		bool _isParallel;
	};
}