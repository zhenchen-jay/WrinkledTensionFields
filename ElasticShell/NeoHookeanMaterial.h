#pragma once
#include "ElasticShellMaterial.h"


/*
* Neo-Hookean nonlinear material model, with energy density
* W = beta/2.0 (tr[M] - 2 - log det M) + alpha/2.0 (log det M / 2)^2
* where M = gbar^-1 g, and g and gbar are the current and rest metrics of the 
* shell volume (which vary in the thickness direction as defined by the surface 
* fundamental forms).
*/

class NeoHookeanMaterial : public ElasticShellMaterial
{
public:
    virtual double stretchingEnergy(
        const MeshConnectivity &mesh,
        const Eigen::MatrixXd &curPos,
        double lameAlpha, double lameBeta, double thickness,
        const Eigen::Matrix2d &abar,
        int face,
        Eigen::Matrix<double, 1, 9> *derivative, // F(face, i)
        Eigen::Matrix<double, 9, 9> *hessian,
        bool isLocalProj = false) override;

    virtual double bendingEnergy(
        const MeshConnectivity &mesh,
        const Eigen::MatrixXd &curPos,
        const Eigen::VectorXd &edgeDOFs,
        double lameAlpha, double lameBeta, double thickness,
        const Eigen::Matrix2d &abar, const Eigen::Matrix2d &bbar,
        int face,
        const SecondFundamentalFormDiscretization &sff,
        Eigen::MatrixXd *derivative, // F(face, i), then the three vertices opposite F(face,i), then the thetas on oppositeEdge(face,i)
        Eigen::MatrixXd *hessian,
        bool isLocalProj = false) override;
};