#include "SecondFundamentalFormDiscretization.h"
#include <iostream>
#include <random>
#include "../MeshLib/MeshConnectivity.h"

void SecondFundamentalFormDiscretization::testSecondFundamentalForm(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
    double eps = 1e-6;
    MeshConnectivity mesh(F);
    int nfaces = mesh.nFaces();
    int nedges = mesh.nEdges();
    int edgedofs = numExtraDOFs();

    Eigen::VectorXd extraDOFs(edgedofs * nedges);
    extraDOFs.setConstant(0.3);
    int ntests = 100;

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uni(0, nfaces-1);

    std::cout << "Testing " << ntests << " random faces" << std::endl;
    for(int i=0; i<ntests; i++)
    {
        int face = uni(rng);
        Eigen::MatrixXd deriv(4, 18 + 3 * edgedofs);
        std::vector<Eigen::MatrixXd > hess;
        Eigen::Matrix2d b = secondFundamentalForm(mesh, V, extraDOFs, face, &deriv, &hess);

        for(int j=0; j<3; j++)
        {
            for(int k=0; k<3; k++)
            {
                Eigen::MatrixXd Vpert(V);
                Vpert(mesh.faceVertex(face, j), k) += 1e-6;
                Eigen::MatrixXd derivpert(4, 18 + 3 * edgedofs);
                Eigen::Matrix2d bpert = secondFundamentalForm(mesh, Vpert, extraDOFs, face, &derivpert, NULL);
                Eigen::Matrix2d findiff = (bpert-b)/1e-6;
                Eigen::Vector4d findiffvec;
                findiffvec << findiff(0,0), findiff(0,1), findiff(1,0), findiff(1,1);
                Eigen::Vector4d exact = deriv.col(3*j+k);
                std::cout << "q" << j << "[" << k <<"]: " << exact.transpose() << " / " << findiffvec.transpose() << std::endl;

                Eigen::MatrixXd findiffhess = (derivpert-deriv)/1e-6;
                for(int l=0; l<4; l++)
                {
                    std::cout << " hess[" << l << "]: " << hess[l].col(3*j+k).transpose() << " / " << findiffhess.row(l) << std::endl;;
                }
            }
            int edge = mesh.faceEdge(face, j);
            int ofaceidx = 0;
            if(mesh.edgeFace(edge, ofaceidx) == face)
                ofaceidx = 1;
            if(mesh.edgeFace(edge, ofaceidx) == -1)
                continue;
            int pidx = mesh.edgeOppositeVertex(edge, ofaceidx);
            for(int k=0; k<3; k++)
            {
                Eigen::MatrixXd Vpert(V);
                Vpert(pidx, k) += 1e-6;
                Eigen::MatrixXd derivpert(4, 18 + 3 * edgedofs);
                Eigen::Matrix2d bpert = secondFundamentalForm(mesh, Vpert, extraDOFs, face, &derivpert, NULL);
                Eigen::Matrix2d findiff = (bpert-b)/1e-6;
                Eigen::Vector4d findiffvec;
                findiffvec << findiff(0,0), findiff(0,1), findiff(1,0), findiff(1,1);
                Eigen::Vector4d exact = deriv.col(9 + 3*j+k);
                std::cout << "p" << j << "[" << k <<"]: " << exact.transpose() << " / " << findiffvec.transpose() << std::endl;

                Eigen::MatrixXd findiffhess = (derivpert-deriv)/1e-6;
                for(int l=0; l<4; l++)
                {
                    std::cout << " hess[" << l << "]: " << hess[l].col(9 + 3*j+k).transpose() << " / " << findiffhess.row(l) << std::endl;;
                }
            }
            Eigen::VectorXd extrapert(extraDOFs);
            for (int k = 0; k < edgedofs; k++)
            {
                extrapert[edgedofs*edge + k] += 1e-6;
                Eigen::MatrixXd derivpert(4, 18 + 3 * edgedofs);
                Eigen::Matrix2d bpert = secondFundamentalForm(mesh, V, extrapert, face, &derivpert, NULL);
                Eigen::Matrix2d findiff = (bpert - b) / 1e-6;
                Eigen::Vector4d findiffvec;
                findiffvec << findiff(0, 0), findiff(0, 1), findiff(1, 0), findiff(1, 1);
                Eigen::Vector4d exact = deriv.col(18 + edgedofs * j + k);
                std::cout << "theta[" << j << "]: " << exact.transpose() << " / " << findiffvec.transpose() << std::endl;

                Eigen::MatrixXd findiffhess = (derivpert - deriv) / 1e-6;
                for (int l = 0; l < 4; l++)
                {
                    std::cout << " hess[" << l << "]: " << hess[l].col(18 + j).transpose() << " / " << findiffhess.row(l) << std::endl;;
                }
            }
        }
    }
}