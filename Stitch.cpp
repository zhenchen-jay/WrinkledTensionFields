#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>
#include <igl/remove_unreferenced.h>
#include "SecondFundamentalForm/MidedgeAverageFormulation.h"
#include "Stitch.h"
#include "MeshConnectivity.h"
#include <igl/cotmatrix.h>
#include <igl/boundary_facets.h>
#include <igl/boundary_loop.h>

void stitchMeshes(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& stitchedV, Eigen::MatrixXi& stitchedF, std::vector<Eigen::Vector3i> &bnd_edges, Eigen::VectorXi& newIndex)
{
    std::map<int, int> parentvert;
    int nverts = V.rows();
    std::vector<bool> seam(nverts,false);
    std::vector<std::vector<int> > bls;
       
    igl::boundary_loop(F, bls);
       
    std::vector<int> boundVerts;
       
    for(auto &b : bls)
        for(int i = 0; i < b.size(); i++)
        {
            boundVerts.push_back(b[i]);
        }
    
    for(int i = 0; i < boundVerts.size(); i++)
    {
        for(int j = 0; j < boundVerts.size(); j++)
        {
            if(j == i)
                continue;
            int v1 = boundVerts[i];
            int v2 = boundVerts[j];
            if((V.row(v1) - V.row(v2)).norm() < 1e-5)
            {
                seam[v1] = true;
                seam[v2] = true;
                auto it = parentvert.find(v1);
                
                if(it == parentvert.end())
                    parentvert[v2] = v1;
                else
                    parentvert[v2] = it->second;
            }
        }
    }

    int nfaces = F.rows();
    Eigen::MatrixXi mappedF = F;

    for (int i = 0; i < nfaces; i++)
    { 
        for (int j = 0; j < 3; j++)
        {
            auto it = parentvert.find(mappedF(i, j));
            if (it != parentvert.end())
                mappedF(i,j) = it->second;
        }
    }
 

    Eigen::VectorXi I;
    Eigen::VectorXi J;
    igl::remove_unreferenced(V, mappedF, stitchedV, stitchedF, I, J);

    std::cout << "number of verts over all patches :  " << V.rows() << " number of verts after stitching : " << stitchedV.rows() << std::endl;
    std::cout << "number of faces over all patches :  " << F.rows() << " number of faces after stitching : " << stitchedF.rows() << std::endl;
    
    nverts = stitchedV.rows();
   
    newIndex = I;
    for (int i = 0; i < newIndex.size(); i++)
    {
        if (newIndex(i) == -1)
        {
            auto it = parentvert.find(i);
            newIndex(i) = I(it->second);
        }
    }

    //compute boundary
 
    Eigen::MatrixXi bndV;
    Eigen::MatrixXi bndF;
    Eigen::MatrixXi OppBndV;
    igl::boundary_facets(F,bndV, bndF, OppBndV);
    bnd_edges.clear();
    for (int i = 0; i < bndV.rows(); i++)
    {
        if (!(seam[bndV(i,0)] && seam[bndV(i,1)])) // both vertices are not on the seam then it's still a boundary
            bnd_edges.push_back(Eigen::Vector3i(bndV(i,0), bndV(i,1), F(bndF(i),OppBndV(i)))); 
    }

}


void stitchMeshesWithTol(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXd &stitchedV, Eigen::MatrixXi &stitchedF, double tol)
{
    std::map<int, int> parentvert;
    int nverts = V.rows();
    std::vector<bool> seam(nverts,false);
    std::vector<std::vector<int> > bls;
       
    igl::boundary_loop(F, bls);
       
    std::vector<int> boundVerts;
       
    for(auto &b : bls)
        for(int i = 0; i < b.size(); i++)
        {
            boundVerts.push_back(b[i]);
        }
       
    for(int i = 0; i < boundVerts.size(); i++)
    {
        for(int j = 0; j < boundVerts.size(); j++)
        {
            if(j == i)
                continue;
            int v1 = boundVerts[i];
            int v2 = boundVerts[j];
            if((V.row(v1) - V.row(v2)).norm() < 1e-5)
            {
                seam[v1] = true;
                seam[v2] = true;
                auto it = parentvert.find(v1);
                
                if(it == parentvert.end())
                    parentvert[v2] = v1;
                else
                    parentvert[v2] = it->second;
            }
        }
    }

    int nfaces = F.rows();
    Eigen::MatrixXi mappedF = F;

    for (int i = 0; i < nfaces; i++)
    { 
        for (int j = 0; j < 3; j++)
        {
            auto it = parentvert.find(mappedF(i, j));
            if (it != parentvert.end())
                mappedF(i,j) = it->second;
        }
    }
 

    Eigen::VectorXi I;
    Eigen::VectorXi J;
    igl::remove_unreferenced(V, mappedF, stitchedV, stitchedF, I, J);

    std::cout << "number of verts over all patches :  " << V.rows() << " number of verts after stitching : " << stitchedV.rows() << std::endl;
    std::cout << "number of faces over all patches :  " << F.rows() << " number of faces after stitching : " << stitchedF.rows() << std::endl;
}

void testStitchMeshes()
{
   /* Eigen::MatrixXd V1(8,3);
    Eigen::MatrixXi F1(4,3);
    Eigen::MatrixXd nV1;
    Eigen::MatrixXi nF1;
    V1.setZero();
    V1.row(0) << 0, 0, 0;
    V1.row(1) << 1, 0, 0;
    V1.row(2) << 0, 1, 0;
    V1.row(3) << 1, 1, 0;
    Eigen::RowVector3d shift;
    shift << -1, 0 ,0;
    for (int j = 0; j < 4; j++)
        V1.row(j+4) = V1.row(j) + shift;
    F1 << 0, 1, 3, 0, 3, 2, 4, 5, 7, 4, 7, 6;
    SimulationSetup setup1;
    MidedgeAverageFormulation sff;
    MeshConnectivity mesh1(F1);
    setup1.buildRestFundamentalForms(mesh1, V1, sff);
    stitchMeshes(V1,F1,nV1,nF1,setup1);

    Eigen::MatrixXd nV2;
    Eigen::MatrixXi nF2;
    SimulationSetup setup2;
    MeshConnectivity mesh2(nF1);
    MidedgeAverageFormulation sff2;
    setup2.buildRestFundamentalForms(mesh2, nV1, sff2);
    stitchMeshes(nV1,nF1,nV2,nF2,setup2);

    std::cout << "laplacian difference : " << (setup1.laplacian - setup2.laplacian).norm() << std::endl;   
    std::cout << "laplace position : " << (setup1.laplacian * nV1).norm() << std::endl;*/

}
