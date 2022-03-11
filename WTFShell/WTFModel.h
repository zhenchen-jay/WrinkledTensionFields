#pragma once
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "../external/cppoptlib/problem.h"
#include "../external/cppoptlib/meta.h"
#include "../SecondFundamentalForm/SecondFundamentalFormDiscretization.h"
#include "../SecondFundamentalForm/MidedgeAverageFormulation.h"
#include "../SecondFundamentalForm/MidedgeAngleSinFormulation.h"
#include "../SecondFundamentalForm/MidedgeAngleTanFormulation.h"
#include "../CommonFunctions.h"

#include "WTFSetup.h"
#include "WTFState.h"

using namespace WTF;
using namespace cppoptlib;
using Eigen::VectorXd;

class WTFModel : public Problem<double>
{
  public:
	using typename cppoptlib::Problem<double>::Scalar;
	using typename cppoptlib::Problem<double>::TVector;
	using typename cppoptlib::Problem<double>::THessian;
	
	void initialization(const WTFSetup setup, const WTFState state, std::map<int, double> *clampedDOFs, std::string filePath, bool isUsePosHess = false, bool isParallel = false);
	void convertCurState2Variables(const WTFState curState, TVector &x);
	void convertVariables2CurState(const TVector x, WTFState &curState);

	double value(const TVector &x);
	double stretchingValue(const TVector &x);
	double bendingValue(const TVector &x);

	void gradient(const TVector &x, TVector &grad);
	void hessian(const TVector& x, THessian& hessian);     
	
	void save(int curIterations, TimeCost curTimeCost, double stepSize, double oldEnergy, double curEnergy, double fDelta, double ampDelta, double dphiDelta, double reg, double* projGradNorm, bool PSDHess, bool isSaveAllIntermediate = true);

	Eigen::VectorXd getFullDir(const TVector &dir) { return _projM.transpose() * dir; }
	int nFreeAmp() {return _freeAmp; }
	int nFreePhi() {return _freePhi; }
	int nFreeDphi() { return _freeDphi; }

	Eigen::SparseMatrix<double> buildIntegrabilityConstraints();

	void testValueAndGradient(const TVector &x);
	void testGradientAndHessian(const TVector& x);

	void setProjM(std::map<int, double> *clampedDOFs);
	Eigen::VectorXd fullGradient(const TVector &grad) {return _projM.transpose() * grad;}

	bool isC2() { return true; }

	double getMaxStep(const TVector& x, const TVector& dir, double step) 
	{
		return 1e15;
	}
	

	// testing whether near the optimality, the constrainted hessian is near PD
	void checkPD4ConstraintedHess(const TVector &x);

	Eigen::VectorXd getProjectedGradient(const TVector &x);

private:
	void getUpsampledWrinkledMesh(Eigen::MatrixXd& NV, Eigen::MatrixXi& NF, Eigen::MatrixXd &upsampledTFTV, Eigen::MatrixXi &upsampledTFTF, Eigen::MatrixXd &soupPhiV, Eigen::MatrixXi &soupPhiF, Eigen::MatrixXd &soupProblemV, Eigen::MatrixXi &soupProblemF, Eigen::VectorXd &upsampledAmp, Eigen::VectorXd &soupPhi);

public:
	WTFSetup _setup;
	WTFState _state;
	//Eigen::MatrixXd _prevPos;
	Eigen::SparseMatrix<double> _projM;
	Eigen::MatrixXd _reducedBasis;
	std::map<int, double> _clampedDOFs;

	bool _isUseD2phi;
	bool _isUsePosHess;

	bool _isParallel;

	int _freeAmp;
	int _freePhi;
	int _freeDphi;

	double _reg;

	std::string _filePrefix;
};
