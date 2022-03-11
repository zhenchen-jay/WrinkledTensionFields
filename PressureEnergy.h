#ifndef PRESSUREENERGY_H
#define PRESSUREENERGY_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <vector>

double pressureEnergy(
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXd &curPos,
    double pressure,
    Eigen::VectorXd *derivative,
    std::vector<Eigen::Triplet<double> > *hessian,
    Eigen::VectorXd center = Eigen::Vector3d::Zero(),
    bool isProjHess = false, 
    bool isParallel = false);

double pressureEnergyPerface(const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& curPos,
    double pressure,
    int face,
    Eigen::VectorXd* derivative,
    Eigen::MatrixXd* hessian,
    Eigen::VectorXd center = Eigen::Vector3d::Zero(),
    bool isProjHess = false);

void testPressureEnergy(
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& curPos,
    double pressure,
    Eigen::VectorXd center = Eigen::Vector3d::Zero()
);

void testPressureEnergyPerface(
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& curPos,
    double pressure,
    int face,
    Eigen::VectorXd center = Eigen::Vector3d::Zero()
);

#endif
