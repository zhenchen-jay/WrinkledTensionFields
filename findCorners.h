#ifndef FINDCORNERS_H
#define FINDCORNERS_H

#include <Eigen/Core>
#include <set>

void findCorners(const Eigen::MatrixXd& V2D, const Eigen::MatrixXi& F2D, const Eigen::MatrixXd &V3D, const Eigen::MatrixXi& F3D, std::set<int> &corners);

#endif