#ifndef PENALTYENERGY_H
#define PENALTYENERGY_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include "../ElasticShell/ElasticSetup.h"
#include "../Obstacle.h"
#include "../CommonFunctions.h"

double penaltyEnergy_Sphere(
    const Eigen::MatrixX3d &V, // Current (deformed) positions of the vertices
    double penaltyK, // penalty stiffness
    Eigen::Vector3d ctr, // radius of the ball centered at (0,0,-r)
    double r, // radius of the ball centered at (0,0,-r)
    Eigen::VectorXd *dEnergy, // If not null, the energy gradient will be written to this vector
    std::vector<Eigen::Triplet<double> > *hEnergy, // If not null, the energy Hessian will be written to this vector (in sparse matrix form)
    bool isProjHess = false
);

double springEnergy_Sphere(
    const Eigen::MatrixX3d &V, // Current (deformed) positions of the vertices
    const std::set<int> &collidingPoints, //points that collides with the sphere at least once
    double penaltyK, // penalty stiffness
    Eigen::Vector3d ctr, // radius of the ball centered at (0,0,-r)
    double r, // radius of the ball centered at (0,0,-r)
    Eigen::VectorXd *dEnergy, // If not null, the energy gradient will be written to this vector
    std::vector<Eigen::Triplet<double> > *hEnergy, // If not null, the energy Hessian will be written to this vector (in sparse matrix form)}
    bool isProjHess = false
);

double penaltyEnergy_Obstacle(
    const Eigen::MatrixXd &V, // Current (deformed) positions of the vertices
    const Eigen::MatrixXi &F, // Current (deformed) positions of the vertices
    const Eigen::MatrixXd &obs_V,
    const Eigen::MatrixXi &obs_F,
    double penaltyK, // penalty stiffness
    Eigen::VectorXd *dE, // If not null, the energy gradient will be written to this vector
    std::vector<Eigen::Triplet<double> > *hE,
    bool isProjHess = false);
//    Eigen::MatrixXd &LocalNearestP,
//    Eigen::VectorXd &LocalMinDis);

double penaltyEnergy_prescribed(
    const Eigen::MatrixXd &V, // Current (deformed) positions of the vertices
    const Eigen::MatrixXi &F,
    double penaltyK, // penalty stiffness
    Eigen::MatrixXd LocalNearestP,
    Eigen::VectorXd LocalMinDis);


double penaltyForce_VertexFace(
    const Eigen::MatrixXd &cloth_start,
    const ElasticSetup &setup,
    Eigen::VectorXd *dE, // If not null, the energy gradient will be written to this vector
    std::vector<Eigen::Triplet<double> > *hE, // If not null, the energy Hessian will be written to this vector (in sparse matrix form),
    bool isProjHess);

double vertexFaceEnergy(
    const Eigen::Vector3d& p,
    const Eigen::Vector3d& q0,
    const Eigen::Vector3d& q1,
    const Eigen::Vector3d& q2,
    const double &threshold,
    const double &coeff,
    Eigen::Vector3d *deriv,
    Eigen::Matrix3d *hess
);

void testPenaltyGradient(
    const Eigen::MatrixXd &V,
    const ElasticSetup &setup);

void testVertexFaceEnergy();

#endif
