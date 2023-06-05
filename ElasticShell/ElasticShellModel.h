#pragma once 
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "ElasticSetup.h"
#include "ElasticState.h"
#include "../SecondFundamentalForm/SecondFundamentalFormDiscretization.h"
#include "../CommonFunctions.h"

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

class ElasticShellModel
{
public:
    bool initialization(const ElasticSetup setup, const ElasticState initialGuess, std::string filePrefix, bool posHess = true, bool isParallel = true);
    void convertCurState2Variables(const ElasticState curState, Eigen::VectorXd& x);
    void convertVariables2CurState(const Eigen::VectorXd x, ElasticState& curState);

    double value(const Eigen::VectorXd& x);
    double stretchingValue(const Eigen::VectorXd& x);
    double bendingValue(const Eigen::VectorXd& x);
    double penaltyValue(const Eigen::VectorXd& x);

    void gradient(const Eigen::VectorXd& x, Eigen::VectorXd& grad);
    void hessian(const Eigen::VectorXd& x, Eigen::SparseMatrix<double>& hessian);

    //max step before touching the obstacles
    double getMaxStep(const Eigen::VectorXd& x, const Eigen::VectorXd& dir, double step);
    void testMaxStep();

    void save(int curIterations, TimeCost curTimeCost, double stepsize, double oldEnergy, double curEnergy, double gradnorm, double dirnorm, double reg, bool PSDHess);

    Eigen::VectorXd getFullDir(const Eigen::VectorXd& dir)
    {
        Eigen::VectorXd ret;
        _proj.unprojectVector(dir, ret);
        return ret;
    }

    Eigen::SparseMatrix<double> buildLinearConstraints();

    void testValueAndGradient(const Eigen::VectorXd& x);
    void testGradientAndHessian(const Eigen::VectorXd& x);

    void setProjM();
    Eigen::VectorXd fullGradient(const Eigen::VectorXd& grad)
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
//    double _interp;
    bool _isC2;
    bool _isUsePosHess;
    bool _isParallel;
};
