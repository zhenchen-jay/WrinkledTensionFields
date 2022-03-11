#pragma once

#include "ElasticSetup.h"
#include "StVKMaterial.h"
#include "NeoHookeanMaterial.h"
#include "StVKTensionFieldMaterial.h"
#include "ElasticShellMaterial.h"

double elasticStretchingEnergy(
    const MeshConnectivity &mesh,
    const Eigen::MatrixXd &curPos,
    const Eigen::VectorXd &extraDOFs,
    double lameAlpha, double lameBeta, double thickness,
    const std::vector<Eigen::Matrix2d> &abars, 
    const SecondFundamentalFormDiscretization& sff,
    ElasticShellMaterial &mat,
    Eigen::VectorXd *derivative, // positions, then thetas
    std::vector<Eigen::Triplet<double> > *hessian,
    bool isLocalProj,
    bool isParallel = false);


double elasticBendingEnergy(
    const MeshConnectivity &mesh,
    const Eigen::MatrixXd &curPos,
    const Eigen::VectorXd &extraDOFs,
    double lameAlpha, double lameBeta, double thickness,
    const std::vector<Eigen::Matrix2d> &abars, 
    const std::vector<Eigen::Matrix2d> &bbars,
    const SecondFundamentalFormDiscretization &sff,
    ElasticShellMaterial &mat,
    Eigen::VectorXd *derivative, // positions, then thetas
    std::vector<Eigen::Triplet<double> > *hessian,
    bool isLocalProj,
    bool isParallel = false);

double quadraticBendingEnergy(
    const ElasticSetup &setup,
    const Eigen::MatrixXd &curPos,
    double lameAlpha,
    Eigen::VectorXd *derivative, 
    std::vector<Eigen::Triplet<double> > *hessianCoeff);

void testElasticStretchingEnergy(const MeshConnectivity& mesh,
    const Eigen::MatrixXd& curPos,
    const Eigen::VectorXd& extraDOFs,
    double lameAlpha, double lameBeta, double thickness,
    const std::vector<Eigen::Matrix2d>& abars,
    const SecondFundamentalFormDiscretization& sff,
    ElasticShellMaterial& mat,
    bool isParallel = false);

void testElasticBendingEnergy(const MeshConnectivity& mesh,
    const Eigen::MatrixXd& curPos,
    const Eigen::VectorXd& extraDOFs,
    double lameAlpha, double lameBeta, double thickness,
    const std::vector<Eigen::Matrix2d>& abars,
    const std::vector<Eigen::Matrix2d>& bbars,
    const SecondFundamentalFormDiscretization& sff,
    ElasticShellMaterial& mat,
    bool isParallel = false);

void testQuadraticBendingEnergy(const MeshConnectivity& mesh,
    const Eigen::MatrixXd& curPos,
    const Eigen::VectorXd& extraDOFs,
    double lameAlpha, double lameBeta, double thickness,
    const std::vector<Eigen::Matrix2d>& abars,
    const SecondFundamentalFormDiscretization& sff,
    ElasticShellMaterial& mat,
    bool isParallel = false);
