#pragma once
#include "ComplexLoop.h"

class ComplexLoopZuenko : public ComplexLoop	// We modify the Loop.h
{
public:
    void virtual Subdivide(const Eigen::VectorXd& omega, const std::vector<std::complex<double>>& zvals, Eigen::VectorXd& omegaNew, std::vector<std::complex<double>>& upZvals, int level) override;

private:
    std::complex<double> interpZ(const std::vector<std::complex<double>>& zList, const std::vector<Eigen::Vector3d>& gradThetaList, std::vector<double>& coords, const std::vector<Eigen::Vector3d>& pList);

    std::complex<double> computeZandGradZ(const Eigen::VectorXd& omega, const std::vector<std::complex<double>>& zvals, int fid, const Eigen::Vector3d& bary, Eigen::Vector3cd* gradz);
    Eigen::Vector3d computeGradThetaFromOmegaPerface(const Eigen::VectorXd& omega, int fid, int vInF);
    Eigen::Vector3d computeBaryGradThetaFromOmegaPerface(const Eigen::VectorXd& omega, int fid, const Eigen::Vector3d& bary);

    void updateLoopedZvals(const Eigen::VectorXd& omega, const std::vector<std::complex<double>>& zvals, std::vector<std::complex<double>>& upZvals);
};
