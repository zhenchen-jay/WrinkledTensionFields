#include <igl/boundary_facets.h>
#include <igl/boundary_loop.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

#include "../MeshLib/GeometryDerivatives.h"
#include "../MeshLib/MeshConnectivity.h"
#include "../CommonFunctions.h"

#include "ElasticSetup.h"

void ElasticSetup::buildRestFundamentalForms()
{
    MeshConnectivity restMesh(restF);
    int nfaces = restMesh.nFaces();
    abars.resize(nfaces);
    bbars.resize(nfaces);

    sff->initializeExtraDOFs(restEdgeDOFs, restMesh, restV);

    for (int i = 0; i < nfaces; i++)
    {
        abars[i] = firstFundamentalForm(restMesh, restV, i, NULL, NULL);
        if (restFlat)
            bbars[i].setZero();
        else
            bbars[i] = sff->secondFundamentalForm(restMesh, restV, restEdgeDOFs, i, NULL, NULL);
    }
}

void ElasticSetup::computeVertArea(const int nverts, const Eigen::MatrixXi& stitchedF)
{
    vertArea.resize(nverts);
    std::fill(vertArea.begin(), vertArea.end(), 0);
    for (int i = 0; i < stitchedF.rows(); i++)
    {
        double area = 0.5 * sqrt(abars[i].determinant()) / 3.0;
        for (int j = 0; j < 3; j++)
            vertArea[stitchedF(i, j)] += area;
    }
}


void ElasticSetup::computeLaplacian(const Eigen::MatrixXd restV, const Eigen::MatrixXi restF, const std::vector<Eigen::Vector3i>& bnd_edges, const Eigen::VectorXi& newIndex, const double nverts)
{
    laplacian.resize(nverts, nverts);
    Eigen::SparseMatrix<double> restLaplacian;
    igl::cotmatrix(restV, restF, restLaplacian);
    for (int i = 0; i < bnd_edges.size(); i++)
    {
        int vidx0 = bnd_edges[i](0);
        int vidx1 = bnd_edges[i](1);
        int vidx2 = bnd_edges[i](2);

        Eigen::Vector3d v0 = restV.row(vidx0); // edge vertex
        Eigen::Vector3d v1 = restV.row(vidx1); // edge vertex
        Eigen::Vector3d v2 = restV.row(vidx2); // edge opposite vertex
        double cotan12 = cotan(v0, v1, v2);
        double cotan02 = cotan(v1, v0, v2);
        restLaplacian.coeffRef(vidx0, vidx0) += 0.5 * cotan02;
        restLaplacian.coeffRef(vidx0, vidx1) += 0.5 * cotan12;
        restLaplacian.coeffRef(vidx0, vidx2) -= 0.5 * (cotan02 + cotan12);
        restLaplacian.coeffRef(vidx1, vidx0) += 0.5 * cotan02;
        restLaplacian.coeffRef(vidx1, vidx1) += 0.5 * cotan12;
        restLaplacian.coeffRef(vidx1, vidx2) -= 0.5 * (cotan02 + cotan12);
    }


    std::vector<Eigen::Triplet<double> > lapList;
    for (int k = 0; k < restLaplacian.outerSize(); k++)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(restLaplacian, k); it; ++it)
        {
            lapList.push_back(Eigen::Triplet<double>(newIndex(it.row()), newIndex(it.col()), it.value()));
        }
    }

    laplacian.setFromTriplets(lapList.begin(), lapList.end());
}
