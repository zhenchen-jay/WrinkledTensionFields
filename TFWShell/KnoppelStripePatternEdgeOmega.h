#pragma once
#include "../MeshLib/MeshConnectivity.h"
#include "../CommonFunctions.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

void computeEdgeMatrix(const MeshConnectivity &mesh, const Eigen::VectorXd& edgeW, const Eigen::VectorXd& edgeWeight, const int nverts, Eigen::SparseMatrix<double>& A);

void roundZvalsFromEdgeOmega(const MeshConnectivity &mesh, const Eigen::VectorXd& edgeW, const Eigen::VectorXd& edgeWeight, const Eigen::VectorXd& vertArea, const int nverts, std::vector<std::complex<double>>& zvals);