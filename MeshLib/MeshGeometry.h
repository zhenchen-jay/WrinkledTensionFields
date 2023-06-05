#pragma once
#include <Eigen/Core>
#include <vector>
#include "MeshConnectivity.h"

class MeshGeometry
{
public:
    MeshGeometry();
    MeshGeometry(const Eigen::MatrixXd &V, MeshConnectivity &mesh);

public:
    Eigen::MatrixXd cDiffs; // cDiffs.row(2*i) is the vector in barycentric coordinates of face edgeFace(i,0) from centroid on face edgeFace(i,0) to that on face edgeFace(i,1) after rotating triangle edgeFace(i,1) to edgeFace(i,0)'s tangent plane
                            // cDiffs.row(2*i+1) is the same but with the roles of the faces reversed
    std::vector<Eigen::Matrix<double, 3, 2> > Bs; // Basis vectors in ambient coordinates for the barycentric coordinates on faces
    Eigen::MatrixXd Ts;     // Transition matrices. Ts.block<2,2>(2*i,0) maps vectors from barycentric coordinates of face edgeFace(i,0) to barycentric coordinates of edgeFace(i,1). Ts.block<2,2>(2*i, 2) is opposite.
    Eigen::MatrixXd Js;     // Js.block<2,2>(2*i,0) rotates vectors on face i (in face i's barycentric coordinates) to the perpendicular vector (as measured in ambient space)
    Eigen::MatrixXd faceNormals; // |F| x 3 matrix of face normals
    double averageEdgeLength; // exactly what it says on the tin
};
