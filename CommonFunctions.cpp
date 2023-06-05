#include <Eigen/SPQRSupport>
#include <igl/per_vertex_normals.h>
#include <iostream>
#include <filesystem>
#include "CommonFunctions.h"
#include "MeshLib/MeshGeometry.h"
#include "MeshLib/IntrinsicGeometry.h"
#include "MeshLib/MeshConnectivity.h"

#ifdef __APPLE__
namespace fs = std::__fs::filesystem;
#else
namespace fs = std::filesystem;
#endif

Eigen::MatrixXd lowRankApprox(Eigen::MatrixXd A)
{
    Eigen::MatrixXd posHess = A;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(posHess);
    Eigen::VectorXd evals = es.eigenvalues();

    for (int i = 0; i < evals.size(); i++)
    {
        if (evals(i) < 0)
            evals(i) = 0;
    }
    Eigen::MatrixXd D = evals.asDiagonal();
    Eigen::MatrixXd V = es.eigenvectors();
    posHess = V * D * V.inverse();

    return posHess;
}

Eigen::MatrixXd nullspaceExtraction(const Eigen::SparseMatrix<double> A)
{
	Eigen::SPQR<Eigen::SparseMatrix<double>> solver_SPQR(A.transpose());
	int r = solver_SPQR.rank();

	Eigen::MatrixXd N(A.cols(), A.cols() - r);
	for (int i = 0; i < A.cols() - r; i++)
	{
		Eigen::VectorXd ei(A.cols());
		ei.setZero();
		ei(i + r) = 1;
		N.col(i) = solver_SPQR.matrixQ() * ei;
	}
	std::cout << "rank of A = " << r << std::endl;
	return N;
}

double cotan(const Eigen::Vector3d v0, const Eigen::Vector3d v1, const Eigen::Vector3d v2)
{
    double e0 = (v2 - v1).norm();
    double e1 = (v2 - v0).norm();
    double e2 = (v0 - v1).norm();
    double angle0 = acos((e1 * e1 + e2 * e2 - e0 * e0) / (2 * e1 * e2));
    double cot = 1.0 / tan(angle0);

    return cot;
}


void locatePotentialPureTensionFaces(const std::vector<Eigen::Matrix2d>& abars, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, std::set<int>& potentialPureTensionFaces)
{
	MeshConnectivity mesh(F);
	MeshGeometry curGeo(V, mesh);
	IntrinsicGeometry restGeo(mesh, abars);
	assert(abars.size() == F.rows());
	int nfaces = F.rows();
	int nverts = V.rows();

	Eigen::VectorXd smallestEvals(nfaces);

	std::vector<bool> isPureFaces(nfaces, false);
	std::vector<bool> isPureVerts(nverts, false);

	for (int i = 0; i < nfaces; i++)
	{
		Eigen::Matrix2d abar = abars[i];
		Eigen::Matrix2d a = curGeo.Bs[i].transpose() * curGeo.Bs[i];
		Eigen::Matrix2d diff = a - abar;
		Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix2d> solver(diff, abar);

		smallestEvals(i) = std::min(solver.eigenvalues()[0], solver.eigenvalues()[1]);
	}


	for (int i = 0; i < nfaces; i++)
	{
		if (smallestEvals(i) >= 0)
		{
			isPureFaces[i] = true;
		}
	}

	for (int i = 0; i < nfaces; i++)
	{
		if (isPureFaces[i])
		{
			for (int j = 0; j < 3; j++)
			{
				isPureVerts[mesh.faceVertex(i, j)] = true;
			}
		}
	}

	//// if a face is not pure tension, but all its vertices are on the other pure faces, we think it is pure tension
	//// this is used to fill the mixed state "holes" in the pure tension region
	for (int i = 0; i < nfaces; i++)
	{
		if (!isPureFaces[i])
		{
			bool isAllPureVerts = true;
			for (int j = 0; j < 3; j++)
			{
				int vid = mesh.faceVertex(i, j);
				if (isPureVerts[vid] == false)
					isAllPureVerts = false;
			}
			isPureFaces[i] = isAllPureVerts;
		}
	}

	for (int i = 0; i < nfaces; i++)
	{
		if (isPureFaces[i])
			potentialPureTensionFaces.insert(i);
	}
	// set this to empty means we never relax the integrability constriant.
	// potentialPureTensionFaces.clear();
}

void getPureTensionVertsEdges(const std::set<int>& potentialPureTensionFaces, const Eigen::MatrixXi& F, std::set<int>* pureTensionEdges, std::set<int>* pureTensionVerts)
{
	MeshConnectivity mesh(F);

	std::set<int> tmpPureTensionEdges, tmpPureTensionVerts;
	
	std::set<int> mixedStateEdges;
	// locate the mixed edges: 
	for (int i = 0; i < mesh.nEdges(); i++)
	{
		Eigen::Vector2i faceFlags;
		faceFlags << 0, 0;
		for (int j = 0; j < 2; j++)
		{
			int fid = mesh.edgeFace(i, j);
			if (potentialPureTensionFaces.find(fid) != potentialPureTensionFaces.end())
			{
				faceFlags(j) = 1;
			}
		}
		if (faceFlags(0) + faceFlags(1) == 1)
		{
			if (mixedStateEdges.find(i) == mixedStateEdges.end())
			{
				mixedStateEdges.insert(i);
			}
		}

		else if (faceFlags(0) + faceFlags(1) == 2)
		{
			if (tmpPureTensionEdges.find(i) == tmpPureTensionEdges.end())
			{
				tmpPureTensionEdges.insert(i);
			}
		}
	}
	// mixedStateEdges.clear();
	// tmpPureTensionEdges.clear();
	std::cout << mixedStateEdges.size() << " edges between compressed and pure tension faces. " << tmpPureTensionEdges.size() << " edges are inside the pure tension region" << std::endl;

	for (auto& fid : potentialPureTensionFaces)
	{
		for (int j = 0; j < 3; j++)
		{
			int vid = mesh.faceVertex(fid, j);
			if (tmpPureTensionVerts.find(vid) == tmpPureTensionVerts.end())
			{
				tmpPureTensionVerts.insert(vid);
			}
		}
	}
	// tmpPureTensionVerts.clear();
	std::cout << tmpPureTensionVerts.size() << " vertices are inside the pure tension region (including the boundary)" << std::endl;

	if (pureTensionEdges)
	{
		*pureTensionEdges = tmpPureTensionEdges;
	}
	if (pureTensionVerts)
	{
		*pureTensionVerts = tmpPureTensionVerts;
	}
}


void trivialOffset(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd &offsettedV, double d)
{
	Eigen::MatrixXd VN;
	igl::per_vertex_normals(V, F, VN);
	offsettedV = V  + d * VN;
}

void rigidBodyAlignment(const Eigen::MatrixXd& tarPos, const MeshConnectivity& mesh, const Eigen::MatrixXd &pos, Eigen::Matrix3d &R, Eigen::Vector3d &t)
{
	int nverts = tarPos.rows();
    int nfaces = mesh.nFaces();
    Eigen::VectorXd massVec;
    massVec.resize(nverts);
    massVec.setZero();
    
    Eigen::VectorXd areaList;
    igl::doublearea(tarPos, mesh.faces(), areaList);
    areaList = areaList / 2;
    
    for(int i=0; i < nfaces; i++)
    {
        double faceArea = areaList(i);
        for(int j=0; j<3; j++)
        {
            int vertIdx = mesh.faceVertex(i, j);
            massVec(vertIdx) += faceArea / 3;
        }
    }
    
    massVec = massVec / 3;
    massVec = massVec / massVec.maxCoeff();
    
    Eigen::MatrixXd W = Eigen::MatrixXd::Zero(nverts, nverts);
    W.diagonal() = massVec;
    
    Eigen::Vector3d avePos, aveTarPos;
    avePos.setZero();
    aveTarPos.setZero();
    
    for(int i = 0; i < nverts; i++)
    {
        avePos += massVec(i) * pos.row(i);
        aveTarPos += massVec(i) * tarPos.row(i);
    }
    
    avePos = avePos / massVec.sum();
    aveTarPos = aveTarPos / massVec.sum();
    
    Eigen::MatrixXd onesMat(nverts,1);
    onesMat.setOnes();
    
    Eigen::MatrixXd shiftedPos = pos - onesMat * avePos.transpose();
    Eigen::MatrixXd shiftedTarPos = tarPos - onesMat * aveTarPos.transpose();
    
    Eigen::MatrixXd S = shiftedPos.transpose() * W * shiftedTarPos;
    Eigen::VectorXd b = Eigen::VectorXd::Zero(3);
    
    
    Eigen::BDCSVD<Eigen::MatrixXd> solver(S);
    solver.compute(S, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd U = solver.matrixU();
    Eigen::MatrixXd V = solver.matrixV();
    
    Eigen::MatrixXd middleMat(S.rows(), S.cols());
    middleMat.setIdentity();
    middleMat(S.rows()-1,S.cols()-1) = (V*U.transpose()).determinant();
    
    R = V * middleMat * U.transpose();
    t = aveTarPos - R*avePos;
}

void matToVec(const Eigen::MatrixXd& mat, Eigen::VectorXd& vec)
{
	int nverts = mat.rows();
	vec.resize(3 * nverts);
	for (int i = 0; i < nverts; i++)
	{
		for (int j = 0; j < 3; j++)
			vec(3 * i + j) = mat(i, j);
	}
}

void vecToMat(const Eigen::VectorXd& vec, Eigen::MatrixXd& mat)
{
	int nverts = vec.size() / 3;
	mat.resize(nverts, 3);
	for (int i = 0; i < nverts; i++)
	{
		for (int j = 0; j < 3; j++)
			mat(i, j) = vec(3 * i + j);
	}
}

bool mkdir(const std::string& foldername)
{
	if (!fs::exists(foldername))
	{
		std::cout << "create directory: " << foldername << std::endl;
		if (!fs::create_directory(foldername))
		{
			std::cerr << "create folder failed." << foldername << std::endl;
			return false;
		}
	}
	return true;
}