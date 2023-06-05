#include "ElasticEnergy.h"
#include <tbb/tbb.h>
#include <iostream>

double elasticStretchingEnergy(
    const MeshConnectivity &mesh,
    const Eigen::MatrixXd &curPos,
    const Eigen::VectorXd &extraDOFs,
    double lameAlpha, double lameBeta, double thickness,
    const std::vector<Eigen::Matrix2d> &abars, 
    const SecondFundamentalFormDiscretization& sff,
    ElasticShellMaterial &mat,
    Eigen::VectorXd *derivative, // positions, then thetas
    std::vector<Eigen::Triplet<double> > *hessian,
    bool isLocalProj,
    bool isParallel)
{
    int nfaces = mesh.nFaces();
    int nedges = mesh.nEdges();
    int nverts = curPos.rows();
    int nedgedofs = sff.numExtraDOFs();

    if (derivative)
    {
        derivative->resize(3 * nverts + nedgedofs * nedges);
        derivative->setZero();
    }
    if (hessian)
    {
        hessian->clear();
    }

    double result = 0;
    
    // stretching terms
    auto energies = std::vector<double>(nfaces);
    auto derivs = std::vector<Eigen::Matrix<double, 1, 9>>(nfaces);
    auto hesses = std::vector<Eigen::Matrix<double, 9, 9>>(nfaces);

    if (isParallel)
    {
        auto computeStretching = [&](const tbb::blocked_range<uint32_t>& range)
        {

            for (uint32_t i = range.begin(); i < range.end(); ++i)
            {
                energies[i] = mat.stretchingEnergy(mesh, curPos, lameAlpha, lameBeta, thickness, abars[i], i, derivative ? &derivs[i] : NULL, hessian ? &hesses[i] : NULL, isLocalProj);
            }

        };

        tbb::blocked_range<uint32_t> rangex(0u, (uint32_t)nfaces);
        tbb::parallel_for(rangex, computeStretching);
    }
   
    else
    {
        for (int i = 0; i < nfaces; i++)
        {
            energies[i] = mat.stretchingEnergy(mesh, curPos, lameAlpha, lameBeta, thickness, abars[i], i, derivative ? &derivs[i] : NULL, hessian ? &hesses[i] : NULL, isLocalProj);
        }
    }

    for (int i = 0; i < nfaces; i++)
    {
        result += energies[i];

        if (derivative)
        {
            for (int j = 0; j < 3; j++)
                derivative->segment<3>(3 * mesh.faceVertex(i, j)) += derivs[i].segment<3>(3 * j);
        }
        if (hessian)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 3; m++)
                        {
                            hessian->push_back(Eigen::Triplet<double>(3 * mesh.faceVertex(i, j) + l, 3 * mesh.faceVertex(i, k) + m, hesses[i](3 * j + l, 3 * k + m)));
                        }
                    }
                }
            }
        }
    }

    return result;
    
}

double elasticBendingEnergy(
    const MeshConnectivity& mesh,
    const Eigen::MatrixXd& curPos,
    const Eigen::VectorXd& extraDOFs,
    //const std::set<int> &collapsingEdges, 
    double lameAlpha, double lameBeta, double thickness,
    const std::vector<Eigen::Matrix2d>& abars,
    const std::vector<Eigen::Matrix2d>& bbars,
    const SecondFundamentalFormDiscretization& sff,
    ElasticShellMaterial& mat,
    Eigen::VectorXd* derivative, // positions, then thetas
    std::vector<Eigen::Triplet<double> >* hessian,
    bool isLocalProj,
    bool isParallel)
{
    int nfaces = mesh.nFaces();
    int nedges = mesh.nEdges();
    int nverts = curPos.rows();
    int nedgedofs = sff.numExtraDOFs();

    if (derivative)
    {
        derivative->resize(3 * nverts + nedgedofs * nedges);
        derivative->setZero();
    }
    if (hessian)
    {
        hessian->clear();
    }

    Eigen::MatrixXi F = mesh.faces();


    double result = 0;
    // bending terms
    auto energies = std::vector<double>(nfaces);
    auto derivs = std::vector<Eigen::MatrixXd>(nfaces);
    auto hesses = std::vector<Eigen::MatrixXd>(nfaces);

    if (isParallel)
    {
        auto computeBending = [&](const tbb::blocked_range<uint32_t>& range)
        {
            for (uint32_t i = range.begin(); i < range.end(); ++i)
            {
                Eigen::MatrixXd deriv(1, 18 + 3 * nedgedofs);
                Eigen::MatrixXd hess(18 + 3 * nedgedofs, 18 + 3 * nedgedofs);
                energies[i] = mat.bendingEnergy(mesh, curPos, extraDOFs, lameAlpha, lameBeta, thickness, abars[i], bbars[i], i, sff, derivative ? &deriv : NULL, hessian ? &hess : NULL, isLocalProj);
                if (derivative)
                    derivs[i] = deriv;
                if (hessian)
                    hesses[i] = hess;
            }
        };

        tbb::blocked_range<uint32_t> rangex(0u, (uint32_t)nfaces);
        tbb::parallel_for(rangex, computeBending);
    }

    else
    {
        for (int i = 0; i < nfaces; i++)
        {
            Eigen::MatrixXd deriv(1, 18 + 3 * nedgedofs);
            Eigen::MatrixXd hess(18 + 3 * nedgedofs, 18 + 3 * nedgedofs);
            energies[i] = mat.bendingEnergy(mesh, curPos, extraDOFs, lameAlpha, lameBeta, thickness, abars[i], bbars[i], i, sff, derivative ? &deriv : NULL, hessian ? &hess : NULL, isLocalProj);
            if (derivative)
                derivs[i] = deriv;
            if (hessian)
                hesses[i] = hess;
        }
    }


    for (int i = 0; i < nfaces; i++)
    {
        result += energies[i];
        if (derivative)
        {
            for (int j = 0; j < 3; j++)
            {
                derivative->segment<3>(3 * mesh.faceVertex(i, j)).transpose() += derivs[i].block<1, 3>(0, 3 * j);
                int oppidx = mesh.vertexOppositeFaceEdge(i, j);
                if (oppidx != -1)
                    derivative->segment<3>(3 * oppidx).transpose() += derivs[i].block<1, 3>(0, 9 + 3 * j);
                for (int k = 0; k < nedgedofs; k++)
                {
                    (*derivative)[3 * nverts + nedgedofs * mesh.faceEdge(i, j) + k] += derivs[i](0, 18 + nedgedofs * j + k);
                }
            }
        }
        if (hessian)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 3; m++)
                        {
                            hessian->push_back(Eigen::Triplet<double>(3 * mesh.faceVertex(i, j) + l, 3 * mesh.faceVertex(i, k) + m, hesses[i](3 * j + l, 3 * k + m)));
                            int oppidxk = mesh.vertexOppositeFaceEdge(i, k);
                            if (oppidxk != -1)
                                hessian->push_back(Eigen::Triplet<double>(3 * mesh.faceVertex(i, j) + l, 3 * oppidxk + m, hesses[i](3 * j + l, 9 + 3 * k + m)));
                            int oppidxj = mesh.vertexOppositeFaceEdge(i, j);
                            if (oppidxj != -1)
                                hessian->push_back(Eigen::Triplet<double>(3 * oppidxj + l, 3 * mesh.faceVertex(i, k) + m, hesses[i](9 + 3 * j + l, 3 * k + m)));
                            if (oppidxj != -1 && oppidxk != -1)
                                hessian->push_back(Eigen::Triplet<double>(3 * oppidxj + l, 3 * oppidxk + m, hesses[i](9 + 3 * j + l, 9 + 3 * k + m)));
                        }
                        for (int m = 0; m < nedgedofs; m++)
                        {
                            hessian->push_back(Eigen::Triplet<double>(3 * mesh.faceVertex(i, j) + l, 3 * nverts + nedgedofs * mesh.faceEdge(i, k) + m, hesses[i](3 * j + l, 18 + nedgedofs * k + m)));
                            hessian->push_back(Eigen::Triplet<double>(3 * nverts + nedgedofs * mesh.faceEdge(i, k) + m, 3 * mesh.faceVertex(i, j) + l, hesses[i](18 + nedgedofs * k + m, 3 * j + l)));
                            int oppidxj = mesh.vertexOppositeFaceEdge(i, j);
                            if (oppidxj != -1)
                            {
                                hessian->push_back(Eigen::Triplet<double>(3 * oppidxj + l, 3 * nverts + nedgedofs * mesh.faceEdge(i, k) + m, hesses[i](9 + 3 * j + l, 18 + nedgedofs * k + m)));
                                hessian->push_back(Eigen::Triplet<double>(3 * nverts + nedgedofs * mesh.faceEdge(i, k) + m, 3 * oppidxj + l, hesses[i](18 + nedgedofs * k + m, 9 + 3 * j + l)));
                            }
                        }
                    }
                    for (int m = 0; m < nedgedofs; m++)
                    {
                        for (int n = 0; n < nedgedofs; n++)
                        {
                            hessian->push_back(Eigen::Triplet<double>(3 * nverts + nedgedofs * mesh.faceEdge(i, j) + m, 3 * nverts + nedgedofs * mesh.faceEdge(i, k) + n, hesses[i](18 + nedgedofs * j + m, 18 + nedgedofs * k + n)));
                        }
                    }
                }
            }
        }
    }

    return result;
}


double quadraticBendingEnergy(
    const ElasticSetup &setup,
    const Eigen::MatrixXd &curPos,
    double lameAlpha,
    Eigen::VectorXd *derivative, 
    std::vector<Eigen::Triplet<double> > *hessianCoeff)
{
    int nVerts = curPos.rows();
    int nDOFs = 3 * nVerts + setup.restEdgeDOFs.size();
    if (derivative)
    {
        derivative->resize(nDOFs);
        derivative->setZero();
    }
    if (hessianCoeff)
        hessianCoeff->clear();

    Eigen::VectorXd curV(3 * nVerts);
    for (int i = 0; i < nVerts; i++)
        curV.segment<3>(3 * i) = curPos.row(i);

    Eigen::SparseMatrix<double> newL(3 * nVerts, 3 * nVerts), bendE(3 * nVerts, 3 * nVerts);

    //extend L to 3n x 3n so it applies on all three directions
    std::vector<Eigen::Triplet<double>> Llist;
    for (int k = 0; k < setup.laplacian.outerSize(); k++){
        for (Eigen::SparseMatrix<double>::InnerIterator it(setup.laplacian,k); it; ++it){
            for (int i = 0; i < 3; i++)
                Llist.push_back(Eigen::Triplet<double>(3 * it.row() + i, 3 * it.col() + i, it.value()));
        }
    }
    newL.setFromTriplets(Llist.begin(), Llist.end());

    std::vector<Eigen::Triplet<double> > areaInvList; 
    Eigen::SparseMatrix<double> areaInv(3 * nVerts, 3 * nVerts);
    for(int i=0; i<nVerts; i++){
        for(int j=0; j<3; j++){
            areaInvList.push_back(Eigen::Triplet<double>(3*i+j, 3*i+j, 1.0/setup.vertArea[i]));
        }
    }

    areaInv.setFromTriplets(areaInvList.begin(), areaInvList.end());

    bendE = newL.transpose() * areaInv * newL; 
    double bendingK = lameAlpha / 12.0 * setup.thickness * setup.thickness * setup.thickness;
    double energy = 0.5 * bendingK * curV.transpose() * bendE * curV;
 
    if(derivative){
       (*derivative).segment(0, 3 * nVerts) = bendingK * curV.transpose() * bendE;
    }

    Eigen::SparseMatrix<double> hessian (3 * nVerts, 3 * nVerts);
    if(hessianCoeff){
      hessian = bendingK * bendE;
      for (int k = 0; k < hessian.outerSize(); ++k){
          for (Eigen::SparseMatrix<double>::InnerIterator it(hessian,k); it; ++it)
               hessianCoeff->push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
      }
    }

    return energy;
}

void testElasticStretchingEnergy(const MeshConnectivity& mesh,
    const Eigen::MatrixXd& curPos,
    const Eigen::VectorXd& extraDOFs,
    double lameAlpha, double lameBeta, double thickness,
    const std::vector<Eigen::Matrix2d>& abars,
    const SecondFundamentalFormDiscretization& sff,
    ElasticShellMaterial& mat,
    bool isParallel)
{
    int nfaces = mesh.nFaces();
    int nedges = mesh.nEdges();
    int nverts = curPos.rows();
    int nedgedofs = sff.numExtraDOFs();

    int nDOFs = 3 * nverts + nedgedofs * nedges;

    Eigen::VectorXd perturb = Eigen::VectorXd::Random(nDOFs);

    Eigen::VectorXd deriv;
    std::vector<Eigen::Triplet<double> > hess;

    double energy = elasticStretchingEnergy(mesh, curPos, extraDOFs, lameAlpha, lameBeta, thickness, abars, sff, mat, &deriv, &hess, false, isParallel);
    

    Eigen::SparseMatrix<double> H(nDOFs, nDOFs);
    H.setFromTriplets(hess.begin(), hess.end());

    Eigen::MatrixXd V = curPos;
    std::vector<double> energyVec;
    std::vector<Eigen::VectorXd> derivVec;

    std::cout << "energy = " << energy << ", || g || = " << deriv.norm() << "|| H || = " << H.norm() << std::endl;
    
    for (int k = 4; k <= 10; k++)
    {
        double eps = std::pow(0.1, k);
        for (int i = 0; i < nverts; i++)
        {
            V.row(i) = curPos.row(i) + eps * perturb.segment<3>(3 * i).transpose();
        }
        Eigen::VectorXd deriv1;
        Eigen::VectorXd curExtraDOFs = extraDOFs;
        for (int j = 0; j < nedgedofs * nedges; j++)
        {
            curExtraDOFs(j) = extraDOFs(j) + perturb(3 * nverts + j);
        }
        double energy1 = elasticStretchingEnergy(mesh, V, curExtraDOFs, lameAlpha, lameBeta, thickness, abars, sff, mat, &deriv1, NULL, false, isParallel);
        energyVec.push_back(energy1);
        derivVec.push_back(deriv1);
    }
    
    // energy-deriv check
    std::cout << "energy derivative check. " << std::endl;
    for (int k = 4; k <= 10; k++)
    {
        double eps = std::pow(0.1, k);
        std::cout << "eps: " << eps << std::endl;
        std::cout << "perturbed energy: " << energyVec[k - 4] << ", energy: " << energy << std::endl;
        std::cout << "finite difference: " << (energyVec[k - 4] - energy) / eps << ", directional derivative: " << deriv.dot(perturb) << ", error: " << std::abs((energyVec[k - 4] - energy) / eps - deriv.dot(perturb)) << std::endl;
    }

    // deriv-hess check
    std::cout << "derivative hessian check. " << std::endl;
    for (int k = 4; k <= 10; k++)
    {
        double eps = std::pow(0.1, k);
        std::cout << "eps: " << eps << std::endl;
        std::cout << "perturbed derivative norm: " << derivVec[k - 4].norm() << ", derivative norm: " << deriv.norm() << std::endl;
        std::cout << "finite difference: " << ((derivVec[k - 4] - deriv) / eps).norm() << ", directional derivative: " << (H * perturb).norm() << ", error: " <<((derivVec[k - 4] - deriv) / eps - H * perturb).norm() << std::endl;
    }

}

void testElasticBendingEnergy(const MeshConnectivity& mesh,
    const Eigen::MatrixXd& curPos,
    const Eigen::VectorXd& extraDOFs,
    double lameAlpha, double lameBeta, double thickness,
    const std::vector<Eigen::Matrix2d>& abars,
    const std::vector<Eigen::Matrix2d>& bbars,
    const SecondFundamentalFormDiscretization& sff,
    ElasticShellMaterial& mat,
    bool isParallel)
{
    int nfaces = mesh.nFaces();
    int nedges = mesh.nEdges();
    int nverts = curPos.rows();
    int nedgedofs = sff.numExtraDOFs();

    int nDOFs = 3 * nverts + nedgedofs * nedges;

    Eigen::VectorXd perturb = Eigen::VectorXd::Random(nDOFs);

    Eigen::VectorXd deriv;
    std::vector<Eigen::Triplet<double> > hess;

    std::cout << "|| curPos || = " << curPos.norm() << ", # cur faces = " << mesh.nFaces() << std::endl;
    std::cout << "alpha = " << lameAlpha << ", beta = " << lameBeta << std::endl;
    std::cout << "thickness_min = " << thickness << ", thickness_max = " << thickness << std::endl;
    std::cout << "abar[0] = " << std::endl << abars[0] << std::endl;

    double energy = elasticBendingEnergy(mesh, curPos, extraDOFs, lameAlpha, lameBeta, thickness, abars, bbars, sff, mat, &deriv, &hess, false, isParallel);

    Eigen::SparseMatrix<double> H(nDOFs, nDOFs);
    H.setFromTriplets(hess.begin(), hess.end());

    Eigen::MatrixXd V = curPos;
    std::vector<double> energyVec;
    std::vector<Eigen::VectorXd> derivVec;

    std::cout << "energy = " << energy << ", || g || = " << deriv.norm() << "|| H || = " << H.norm() << std::endl;

    for (int k = 4; k <= 10; k++)
    {
        double eps = std::pow(0.1, k);
        for (int i = 0; i < nverts; i++)
        {
            V.row(i) = curPos.row(i) + eps * perturb.segment<3>(3 * i).transpose();
        }
        Eigen::VectorXd deriv1;
        Eigen::VectorXd curExtraDOFs = extraDOFs;
        for (int j = 0; j < nedgedofs * nedges; j++)
        {
            curExtraDOFs(j) = extraDOFs(j) + eps * perturb(3 * nverts + j);
        }
        double energy1 = elasticBendingEnergy(mesh, V, curExtraDOFs, lameAlpha, lameBeta, thickness, abars, bbars, sff, mat, &deriv1, NULL, false, isParallel);
        energyVec.push_back(energy1);
        derivVec.push_back(deriv1);
    }

    // energy-deriv check
    std::cout << "energy derivative check. " << std::endl;
    for (int k = 4; k <= 10; k++)
    {
        double eps = std::pow(0.1, k);
        std::cout << "eps: " << eps << std::endl;
        std::cout << "perturbed energy: " << energyVec[k - 4] << ", energy: " << energy << std::endl;
        std::cout << "finite difference: " << (energyVec[k - 4] - energy) / eps << ", directional derivative: " << deriv.dot(perturb) << ", error: " << std::abs((energyVec[k - 4] - energy) / eps - deriv.dot(perturb)) << std::endl;
    }

    // deriv-hess check
    std::cout << "derivative hessian check. " << std::endl;
    for (int k = 4; k <= 10; k++)
    {
        double eps = std::pow(0.1, k);
        std::cout << "eps: " << eps << std::endl;
        std::cout << "perturbed derivative norm: " << derivVec[k - 4].norm() << ", derivative norm: " << deriv.norm() << std::endl;
        std::cout << "finite difference: " << ((derivVec[k - 4] - deriv) / eps).norm() << ", directional derivative: " << (H * perturb).norm() << ", error: " << ((derivVec[k - 4] - deriv) / eps - H * perturb).norm() << std::endl;
    }
}

void testQuadraticBendingEnergy(const MeshConnectivity& mesh,
    const Eigen::MatrixXd& curPos,
    const Eigen::VectorXd& extraDOFs,
    double lameAlpha, double lameBeta, double thickness,
    const std::vector<Eigen::Matrix2d>& abars,
    const SecondFundamentalFormDiscretization& sff,
    ElasticShellMaterial& mat,
    bool isParallel)
{

}