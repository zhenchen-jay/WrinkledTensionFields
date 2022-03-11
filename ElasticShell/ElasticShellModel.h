#pragma once 
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "../external/cppoptlib/problem.h"
#include "../external/cppoptlib/meta.h"
#include "ElasticSetup.h"
#include "ElasticState.h"
#include "../SecondFundamentalForm/SecondFundamentalFormDiscretization.h"
#include "../CommonFunctions.h"

using namespace cppoptlib;
using Eigen::VectorXd;

class Projection
{
public:
    Projection() : isPorjNeeded(true) {}
    Projection(const std::vector<bool>& keepDOFs);

    int projDOFs() const { return invdofmap.size(); }
    void projectVector(const Eigen::VectorXd& fullVec, Eigen::VectorXd& projVec) const;
    void unprojectVector(const Eigen::VectorXd& projVec, Eigen::VectorXd& fullVec) const;
    void projectMatrix(std::vector<Eigen::Triplet<double> >& mat) const;

private:
    std::vector<int> dofmap;
    std::vector<int> invdofmap;
    bool isPorjNeeded;
};

class ElasticShellModel : public Problem<double>
{
public:
    using typename cppoptlib::Problem<double>::Scalar;
    using typename cppoptlib::Problem<double>::TVector;
    using typename cppoptlib::Problem<double>::THessian;

    bool initialization(const ElasticSetup setup, const ElasticState initialGuess, std::string filePrefix, double interp = 1.0, bool posHess = true, bool isParallel = true);
    void convertCurState2Variables(const ElasticState curState, TVector& x);
    void convertVariables2CurState(const TVector x, ElasticState& curState);

    double value(const TVector& x);
    double stretchingValue(const TVector& x);
    double bendingValue(const TVector& x);
    double penaltyValue(const TVector& x);

    void gradient(const TVector& x, TVector& grad);
    void hessian(const TVector& x, THessian& hessian);

    //max step before touching the obstacles
    double getMaxStep(const TVector& x, const TVector& dir, double step);
    void testMaxStep();

    void save(int curIterations, TimeCost curTimeCost, double stepsize, double oldEnergy, double curEnergy, double gradnorm, double dirnorm, double reg, bool PSDHess);

    Eigen::VectorXd getFullDir(const TVector& dir)
    {
        Eigen::VectorXd ret;
        _proj.unprojectVector(dir, ret);
        return ret;
    }

    Eigen::SparseMatrix<double> buildLinearConstraints();

    void testValueAndGradient(const TVector& x);
    void testGradientAndHessian(const TVector& x);

    void setProjM();
    Eigen::VectorXd fullGradient(const TVector& grad)
    {
        Eigen::VectorXd ret;
        _proj.unprojectVector(grad, ret);
        return ret;
    }

    bool isC2() { return _isC2; }

    void increasePenaltyStiffness() { _setup.penaltyK * 2.0;  }

public:
    ElasticSetup _setup;
    ElasticState _state;
    Projection _proj;
    double _lameAlpha;
    double _lameBeta;
    std::string _filePrefix;
    double _interp;
    bool _isC2;
    bool _isUsePosHess;
    bool _isParallel;
};
