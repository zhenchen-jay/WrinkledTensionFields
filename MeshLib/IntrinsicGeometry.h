#ifndef INTRINSIC_GEOMETRY
#define INTRINSIC_GEOMETRY

#include <Eigen/Core>
#include <vector>

class MeshConnectivity;

/*
 * Builds data structures for geometric information that depends only on the intrinsic geometry (i.e. the triangle rest metrics)
 * and not on vertex information.
 */
class IntrinsicGeometry
{
public:
    IntrinsicGeometry(const MeshConnectivity &mesh, const std::vector<Eigen::Matrix2d> &abars);

    std::vector<Eigen::Matrix2d> abars; // |F| list of first fundamental forms
    Eigen::MatrixXd Ts;                 // Transition matrices. Ts.block<2,2>(2*i,0) maps vectors from barycentric coordinates of face edgeFace(i,0) to barycentric coordinates of edgeFace(i,1). Ts.block<2,2>(2*i, 2) is opposite.
    Eigen::MatrixXd Js;                 // Js.block<2,2>(2*i,0) rotates vectors on face i (in face i's barycentric coordinates) to the perpendicular vector (under the metric abars[i])
};

#endif