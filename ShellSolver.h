#pragma once
#include "GeometryDerivatives.h"

#include "external/cppoptlib/solver/lbfgssolver.h"
#include "external/cppoptlib/solver/newtondescentsolver.h"
#include "external/cppoptlib/solver/nasoqsovler.h"
#include "WTFShell/WTFModel.h"
#include "WTFShell/WTFSetup.h"
#include "WTFShell/WTFState.h"

#include "ElasticShell/ElasticShellModel.h"
#include "SecondFundamentalForm/SecondFundamentalFormDiscretization.h"

#include "CommonFunctions.h"

struct WTFOptimizationParams
{
	int iterations = 1000;
	double fDelta = 0;
	double gradNorm = 1e-6;
	double xDelta = 0;
	bool isParallel = true;
};

struct FullSimOptimizationParams
{
	size_t iterations = 1000;
	double fDelta = 0;
	double gradNorm = 1e-6;
	double xDelta = 0;
	double interp = 1.0;
	bool isProjH = true;
	bool isParallel = true;
};

// Remark: for dphi version, we only allow entirely fixed dphi 

namespace ShellSolver
{
    void WTFSQPSolver(const WTFSetup& setup, WTFState& curState, std::string filePrefix, const WTFOptimizationParams params);

	void fullSimNewtonStaticSolver(const ElasticSetup& setup, ElasticState& curState, std::string filePrefix, const FullSimOptimizationParams params);
};
