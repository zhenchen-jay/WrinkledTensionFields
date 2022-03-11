#include "MeshGeometry.h"
#include "MeshConnectivity.h"
#include <Eigen/Dense>

MeshGeometry::MeshGeometry()
{
    
}

MeshGeometry::MeshGeometry(const Eigen::MatrixXd &V, MeshConnectivity &mesh)
{
    // compute barycentric matrices and Js
    int nfaces = mesh.nFaces();
    Bs.resize(nfaces);
    Js.resize(2 * nfaces, 2);
    faceNormals.resize(nfaces, 3);

    averageEdgeLength = 0;
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d v0 = V.row(mesh.faces()(i, 0));
        Eigen::Vector3d v1 = V.row(mesh.faces()(i, 1));
        Eigen::Vector3d v2 = V.row(mesh.faces()(i, 2));
        Bs[i].col(0) = v1 - v0;
        Bs[i].col(1) = v2 - v0;

        averageEdgeLength += (v1 - v0).norm();
        averageEdgeLength += (v2 - v1).norm();
        averageEdgeLength += (v0 - v2).norm();

        Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
        n /= n.norm();
        faceNormals.row(i) = n.transpose();

        Eigen::Matrix2d BTB = Bs[i].transpose() * Bs[i];
        Eigen::Matrix<double, 3, 2> ncrossB;
        ncrossB.col(0) = n.cross(v1 - v0);
        ncrossB.col(1) = n.cross(v2 - v0);
        Js.block<2, 2>(2 * i, 0) = BTB.inverse() * Bs[i].transpose() * ncrossB;
    }

    averageEdgeLength /= 3.0 * nfaces;

    // compute cDiffs and transition matrices
    int nedges = mesh.nEdges();
    cDiffs.resize(2 * nedges, 2);
    Ts.resize(2 * nedges, 4);
    for (int edgeidx = 0; edgeidx < nedges; edgeidx++)
    {        
        //collect neighboring face indices
        int face1 = mesh.edgeFace(edgeidx, 0);
        int face2 = mesh.edgeFace(edgeidx, 1);

        if(face1 == -1 || face2 == -1)
            continue;

        int v1 = mesh.edgeVertex(edgeidx, 0);
        int v2 = mesh.edgeVertex(edgeidx, 1);

        Eigen::Vector3d n1 = faceNormals.row(face1).transpose();
        Eigen::Vector3d n2 = faceNormals.row(face2).transpose();
        // collect: (1) the midpoint of the common edge, (2) unit vector in direction of common edge,
        // (3) the face normals, (4) the centroid of neighboring faces

        Eigen::Vector3d midpt = 0.5 * (V.row(v1).transpose() + V.row(v2).transpose());
        Eigen::Vector3d commone = V.row(v2).transpose() - V.row(v1).transpose();
        commone /= commone.norm();
        Eigen::Vector3d centroids[2];
        centroids[0].setZero();
        centroids[1].setZero();

        for (int i = 0; i < 3; i++)
        {
            centroids[0] += V.row(mesh.faces()(face1, i)).transpose();
            centroids[1] += V.row(mesh.faces()(face2, i)).transpose();
        }

        centroids[0] /= 3.0;
        centroids[1] /= 3.0;

        //rotate each centroid into the plane of the opposite triangle and compute ci minus c

        Eigen::Vector3d t1 = n1.cross(commone);
        Eigen::Vector3d t2 = n2.cross(commone);
        Eigen::Vector3d diff2 = centroids[1] - midpt;
        double alpha = commone.dot(diff2);
        double beta = t2.dot(diff2);

        Eigen::Matrix2d BTB1 = Bs[face1].transpose() * Bs[face1];
        Eigen::Matrix2d BTB2 = Bs[face2].transpose() * Bs[face2];

        cDiffs.row(2 * edgeidx) = BTB1.inverse() * Bs[face1].transpose() * (midpt + alpha * commone + beta * t1 - centroids[0]);
        Eigen::Vector3d diff1 = centroids[0] - midpt;
        alpha = commone.dot(diff1);
        beta = t1.dot(diff1);
        cDiffs.row(2 * edgeidx + 1) = BTB2.inverse() * Bs[face2].transpose() * (midpt + alpha*commone + beta * t2 - centroids[1]);

        Eigen::Vector3d e1 = V.row(mesh.faces()(face1, 1)).transpose() - V.row(mesh.faces()(face1, 0)).transpose();
        Eigen::Vector3d e2 = V.row(mesh.faces()(face1, 2)).transpose() - V.row(mesh.faces()(face1, 0)).transpose();

        double alpha1 = commone.dot(e1);
        double beta1 = t1.dot(e1);
        Eigen::Vector3d newe1 = alpha1*commone + beta1 * t2;
        Ts.block<2, 1>(2 * edgeidx, 0) = BTB2.inverse() * Bs[face2].transpose() * newe1;

        double alpha2 = commone.dot(e2);
        double beta2 = t1.dot(e2);
        Eigen::Vector3d newe2 = alpha2*commone + beta2*t2;
        Ts.block<2, 1>(2 * edgeidx, 1) = BTB2.inverse() * Bs[face2].transpose() * newe2;

        e1 = V.row(mesh.faces()(face2, 1)).transpose() - V.row(mesh.faces()(face2, 0)).transpose();
        e2 = V.row(mesh.faces()(face2, 2)).transpose() - V.row(mesh.faces()(face2, 0)).transpose();

        alpha1 = commone.dot(e1);
        beta1 = t2.dot(e1);
        newe1 = alpha1 * commone + beta1 * t1;
        Ts.block<2, 1>(2 * edgeidx, 2) = BTB1.inverse() * Bs[face1].transpose() * newe1;

        alpha2 = commone.dot(e2);
        beta2 = t2.dot(e2);
        newe2 = alpha2*commone + beta2*t1;
        Ts.block<2, 1>(2 * edgeidx, 3) = BTB1.inverse() * Bs[face1].transpose() * newe2;
    }
}
