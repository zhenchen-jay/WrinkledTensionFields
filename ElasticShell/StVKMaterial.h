#pragma once
#include "ElasticShellMaterial.h"

/*
 * St. Venant-Kirchhoff linear material model, with energy density
 * W = alpha/2.0 tr(S)^2 + beta tr(S^2),
 * for strain tensor S = gbar^{-1}(g-gbar), where g and gbar are the current 
 * and rest metrics of the shell volume (which vary in the thickness direction
 * as defined by the surface fundamental forms).
 */

class StVKMaterial : public ElasticShellMaterial
{
public:
    virtual double stretchingEnergy(
        const MeshConnectivity& mesh,
        const Eigen::MatrixXd& curPos,
        double lameAlpha, double lameBeta, double thickness,
        const Eigen::Matrix2d& abar,
        int face,
        Eigen::Matrix<double, 1, 9>* derivative, // F(face, i)
        Eigen::Matrix<double, 9, 9>* hessian,
        bool isLocalProj = false) override;

    virtual double bendingEnergy(
        const MeshConnectivity& mesh,
        const Eigen::MatrixXd& curPos,
        const Eigen::VectorXd& edgeDOFs,
        double lameAlpha, double lameBeta, double thickness,
        const Eigen::Matrix2d& abar, const Eigen::Matrix2d& bbar,
        int face,
        const SecondFundamentalFormDiscretization& sff,
        Eigen::MatrixXd* derivative, // F(face, i), then the three vertices opposite F(face,i), then the thetas on oppositeEdge(face,i)
        Eigen::MatrixXd* hessian,
        bool isLocalProj = false) override;
};