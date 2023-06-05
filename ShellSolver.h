#pragma once
#include "MeshLib/GeometryDerivatives.h"

#include "ElasticShell/ElasticShellModel.h"
#include "TFWShell/TFWModel.h"
#include "SecondFundamentalForm/SecondFundamentalFormDiscretization.h"

#include "CommonFunctions.h"


struct FullSimOptimizationParams
{
	size_t iterations = 1000;
	double fDelta = 0;
	double gradNorm = 1e-6;
	double xDelta = 0;
//	double interp = 1.0;
	bool isProjH = true;
	bool isParallel = true;
	bool printLog = true;
};

struct TFWOptimizationParams
{
	int iterations = 1000;
	double fDelta = 0;
	double gradNorm = 1e-6;
	double xDelta = 0;
	bool isParallel = true;
	bool printLog = true;
};

// Remark: for dphi version, we only allow entirely fixed dphi 

namespace ShellSolver
{
	void fullSimNewtonStaticSolver(const ElasticSetup& setup, ElasticState& curState, std::string filePrefix, const FullSimOptimizationParams& params);
	void TFWSQPSolver(const TFWSetup& setup, TFWState& curState, std::string filePrefix, const TFWOptimizationParams& params);
};
