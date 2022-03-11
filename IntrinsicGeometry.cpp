#include "IntrinsicGeometry.h"
#include "MeshConnectivity.h"
#include <Eigen/Dense>

IntrinsicGeometry::IntrinsicGeometry(const MeshConnectivity &mesh, const std::vector<Eigen::Matrix2d> &abars) : abars(abars)
{
    int nedges = mesh.nEdges();
    int nfaces = mesh.nFaces();

    Js.resize(2 * nfaces, 2);
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Matrix2d Jeuclid;
        Jeuclid << 0, -1,
            1, 0;
        Js.block<2, 2>(2 * i, 0) = std::sqrt(abars[i].determinant()) * abars[i].inverse() * Jeuclid;
    }

    Ts.resize(2 * nedges, 4);
    Ts.setZero();
    for (int i = 0; i < nedges; i++)
    {
        int face1 = mesh.edgeFace(i, 0);
        int face2 = mesh.edgeFace(i, 1);
        if (face1 == -1 || face2 == -1)
            continue;

        int vert1 = mesh.edgeVertex(i, 0);
        int vert2 = mesh.edgeVertex(i, 1);

        // write edge vert1->vert2 in barycentric coordinates on each face
        Eigen::Vector2d barys[3] = { {0,0}, {1,0}, {0,1} };
        Eigen::Vector2d face1e(0, 0);
        Eigen::Vector2d face2e(0, 0);
        for (int j = 0; j < 3; j++)
        {
            if (vert1 == mesh.faceVertex(face1, j))
                face1e -= barys[j];
            else if (vert2 == mesh.faceVertex(face1, j))
                face1e += barys[j];
            if (vert1 == mesh.faceVertex(face2, j))
                face2e -= barys[j];
            else if (vert2 == mesh.faceVertex(face2, j))
                face2e += barys[j];
        }

        Eigen::Vector2d face1eperp = Js.block<2, 2>(2 * face1, 0) * face1e;
        Eigen::Vector2d face2eperp = Js.block<2, 2>(2 * face2, 0) * face2e;
        Eigen::Matrix2d face1basis;
        face1basis.col(0) = face1e;
        face1basis.col(1) = face1eperp;
        Eigen::Matrix2d face2basis;
        face2basis.col(0) = face2e;
        face2basis.col(1) = face2eperp;
        Ts.block<2, 2>(2 * i, 0) = face2basis * face1basis.inverse();
        Ts.block<2, 2>(2 * i, 2) = face1basis * face2basis.inverse();
    }
}
