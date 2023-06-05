#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <vector>
#include <set>

#include "../MeshLib/GeometryDerivatives.h"
#include "../SecondFundamentalForm/SecondFundamentalFormDiscretization.h"

class ElasticShellMaterial
{
    
public:
    ElasticShellMaterial() {}
    virtual ~ElasticShellMaterial() = default;

    virtual double stretchingEnergy(
        const MeshConnectivity &mesh,
        const Eigen::MatrixXd &curPos,
        double lameAlpha, double lameBeta, double thickness,
        const Eigen::Matrix2d &abar,
        int face,
        Eigen::Matrix<double, 1, 9> *derivative, // F(face, i)
        Eigen::Matrix<double, 9, 9> *hessian,
        bool isLocalProj = false) = 0;

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
        bool isLocalProj = false) = 0;
    
    double cotan_v0(const Eigen::Vector3d v0, const Eigen::Vector3d v1, const Eigen::Vector3d v2);
};




