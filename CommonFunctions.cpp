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

Eigen::MatrixXd intrinsicHalfEdgeVec2VertexVec(const Eigen::MatrixXd& v, const Eigen::MatrixXd& pos, const MeshConnectivity& mesh)
{
	int nedges = mesh.nEdges();
	int nverts = pos.rows();
	Eigen::MatrixXd vertOmega(nverts, 3);
	vertOmega.setZero();

	Eigen::MatrixXd vertNormals;
	igl::per_vertex_normals(pos, mesh.faces(), vertNormals);

	Eigen::SparseMatrix<double> A;
	std::vector<Eigen::Triplet<double>> T;

	Eigen::VectorXd edgeVec(2 * nedges);

	for (int i = 0; i < nedges; i++)
	{
		int vid0 = mesh.edgeVertex(i, 0);
		int vid1 = mesh.edgeVertex(i, 1);

		Eigen::Vector3d e = pos.row(vid1) - pos.row(vid0);
		edgeVec.segment<2>(2 * i) = v.row(i);
		for (int j = 0; j < 3; j++)
		{
			T.push_back({ 2 * i, 3 * vid0 + j, e(j) });
			T.push_back({ 2 * i + 1, 3 * vid1 + j, -e(j) });
		}
	}
	A.resize(2 * nedges, 3 * nverts);
	A.setFromTriplets(T.begin(), T.end());

	Eigen::SparseMatrix<double> AT, AAT;
	AT = A.transpose();
	AAT = AT * A;

	Eigen::VectorXd ATb = AT * edgeVec;

	T.clear();
	for (int k=0; k<AAT.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(AAT,k); it; ++it)
		{
			T.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}

	for(int i = 0; i < nverts; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			T.push_back({3 * i + j, 3 * nverts + i, vertNormals(i, j)});
			T.push_back({3 * nverts + i, 3 * i + j, vertNormals(i, j)});
		}
	}

	A.resize(4 * nverts, 4 * nverts);
	A.setFromTriplets(T.begin(), T.end());

	Eigen::VectorXd rhs(4 * nverts);
	rhs.setZero();
	rhs.segment(0, 3 * nverts) = ATb;

	//Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver(A);
	/*Eigen::SPQR<Eigen::SparseMatrix<double>> solver(A);
	Eigen::VectorXd sol = solver.solve(rhs);

	for (int i = 0; i < nverts; i++)
	{
		vertOmega.row(i) = sol.segment<3>(3 * i);
	}*/
	return vertOmega;
}

Eigen::MatrixXd intrinsicHalfEdgeVec2FaceVec(const Eigen::MatrixXd& w, const Eigen::MatrixXd& pos, const MeshConnectivity& mesh)
{
	int nfaces = mesh.nFaces();

	Eigen::MatrixXd faceVec = Eigen::MatrixXd::Zero(nfaces, 3);
	for(int i = 0; i < nfaces; i++)
	{
//        std::cout << "face: " << i  << " of total faces " << nfaces << std::endl;
		for(int j = 0; j < 3; j++)
		{
			int vid = mesh.faceVertex(i, j);

			int eid0 = mesh.faceEdge(i, (j + 1) % 3);
			int eid1 = mesh.faceEdge(i, (j + 2) % 3);

			Eigen::Vector3d e0 = pos.row(mesh.faceVertex(i, (j + 2) % 3)) - pos.row(vid);
			Eigen::Vector3d e1 = pos.row(mesh.faceVertex(i, (j + 1) % 3)) - pos.row(vid);

			int flag0 = 0, flag1 = 0;
			Eigen::Vector2d rhs;

			if (mesh.edgeVertex(eid0, 0) == vid)
			{
				flag0 = 0;
				rhs(0) = w(eid0, 0);
			}
			else
			{
				flag0 = 1;
				rhs(0) = w(eid0, 1);
			}


			if (mesh.edgeVertex(eid1, 0) == vid)
			{
				flag1 = 0;
				rhs(1) = w(eid1, 0);
			}
			else
			{
				flag1 = 1;
				rhs(1) = w(eid1, 1);
			}

			Eigen::Matrix2d I;
			I << e0.dot(e0), e0.dot(e1), e1.dot(e0), e1.dot(e1);
			Eigen::Vector2d sol = I.inverse() * rhs;

			faceVec.row(i) += (sol(0) * e0 + sol(1) * e1) / 3;
		}
	}
	return faceVec;
}

Eigen::MatrixXd intrinsicEdgeVec2FaceVec(const Eigen::VectorXd& w, const Eigen::MatrixXd& pos, const MeshConnectivity& mesh)
{
	int nfaces = mesh.nFaces();

	Eigen::MatrixXd faceVec = Eigen::MatrixXd::Zero(nfaces, 3);
	for (int i = 0; i < nfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int vid = mesh.faceVertex(i, j);

			int eid0 = mesh.faceEdge(i, (j + 1) % 3);
			int eid1 = mesh.faceEdge(i, (j + 2) % 3);

			Eigen::Vector3d e0 = pos.row(mesh.faceVertex(i, (j + 2) % 3)) - pos.row(vid);
			Eigen::Vector3d e1 = pos.row(mesh.faceVertex(i, (j + 1) % 3)) - pos.row(vid);

			int flag0 = 1, flag1 = 1;
			Eigen::Vector2d rhs;

			if (mesh.edgeVertex(eid0, 0) == vid)
			{
				flag0 = 1;
			}
			else
			{
				flag0 = -1;
			}


			if (mesh.edgeVertex(eid1, 0) == vid)
			{
				flag1 = 1;
			}
			else
			{
				flag1 = -1;
			}
			rhs(0) = flag0 * w(eid0);
			rhs(1) = flag1 * w(eid1);

			Eigen::Matrix2d I;
			I << e0.dot(e0), e0.dot(e1), e1.dot(e0), e1.dot(e1);
			Eigen::Vector2d sol = I.inverse() * rhs;

			faceVec.row(i) += (sol(0) * e0 + sol(1) * e1) / 3;
		}
	}
	return faceVec;
}

Eigen::VectorXd faceVec2IntrinsicEdgeVec(const Eigen::MatrixXd& v, const Eigen::MatrixXd& pos, const MeshConnectivity& mesh)
{
	int nedges = mesh.nEdges();
	Eigen::VectorXd edgeOmega(nedges);

	for (int i = 0; i < nedges; i++)
	{
		int fid0 = mesh.edgeFace(i, 0);
		int fid1 = mesh.edgeFace(i, 1);

		int vid0 = mesh.edgeVertex(i, 0);
		int vid1 = mesh.edgeVertex(i, 1);

		Eigen::Vector3d e = pos.row(vid1) - pos.row(vid0);

		if (fid0 == -1)
			edgeOmega(i) = v.row(fid1).dot(e);
		else if (fid1 == -1)
			edgeOmega(i) = v.row(fid0).dot(e);
		else
			edgeOmega(i) = (v.row(fid0) + v.row(fid1)).dot(e) / 2;
	}
	return edgeOmega;
}

Eigen::VectorXd vertexVec2IntrinsicVec(const Eigen::MatrixXd& v, const Eigen::MatrixXd& pos, const MeshConnectivity& mesh)
{
	int nedges = mesh.nEdges();
	Eigen::VectorXd edgeOmega(nedges);

	for (int i = 0; i < nedges; i++)
	{
		int vid0 = mesh.edgeVertex(i, 0);
		int vid1 = mesh.edgeVertex(i, 1);

		Eigen::Vector3d e = pos.row(vid1) - pos.row(vid0);
		edgeOmega(i) = (v.row(vid0) + v.row(vid1)).dot(e) / 2;
	}
	return edgeOmega;
}

Eigen::MatrixXd vertexVec2IntrinsicHalfEdgeVec(const Eigen::MatrixXd& v, const Eigen::MatrixXd& pos, const MeshConnectivity& mesh)
{
	int nedges = mesh.nEdges();
	Eigen::MatrixXd edgeOmega(nedges, 2);

	for (int i = 0; i < nedges; i++)
	{
		int vid0 = mesh.edgeVertex(i, 0);
		int vid1 = mesh.edgeVertex(i, 1);

		Eigen::Vector3d e = pos.row(vid1) - pos.row(vid0);
		edgeOmega(i, 0) = v.row(vid0).dot(e);
		edgeOmega(i, 1) = -v.row(vid1).dot(e);
	}
	return edgeOmega;
}

void computeBaryGradient(const Eigen::Vector3d& P0, const Eigen::Vector3d& P1, const Eigen::Vector3d& P2, const Eigen::Vector3d& bary, Eigen::Matrix3d& baryGrad)
{
	//P = bary(0) * P0 + bary(1) * P1 + bary(2) * P2;
	Eigen::Matrix2d I, Iinv;

	I << (P1 - P0).squaredNorm(), (P1 - P0).dot(P2 - P0), (P2 - P0).dot(P1 - P0), (P2 - P0).squaredNorm();
	Iinv = I.inverse();

	Eigen::Matrix<double, 3, 2> dr;
	dr.col(0) = P1 - P0;
	dr.col(1) = P2 - P0;

	Eigen::Matrix<double, 2, 3> dbary12 = Iinv * dr.transpose();

	baryGrad.row(0) = -dbary12.row(0) - dbary12.row(1);
	baryGrad.row(1) = dbary12.row(0);
	baryGrad.row(2) = dbary12.row(1);
}

Eigen::VectorXd getFaceArea(const Eigen::MatrixXd& V, const MeshConnectivity& mesh)
{
	Eigen::VectorXd faceArea;
	igl::doublearea(V, mesh.faces(), faceArea);
	faceArea /= 2;
	return faceArea;
}

Eigen::VectorXd getEdgeArea(const Eigen::MatrixXd& V, const MeshConnectivity& mesh)
{
	Eigen::VectorXd faceArea = getFaceArea(V, mesh);
	Eigen::VectorXd edgeArea;
	edgeArea.setZero(mesh.nEdges());

	for (int i = 0; i < mesh.nEdges(); i++)
	{
		int f0 = mesh.edgeFace(i, 0);
		int f1 = mesh.edgeFace(i, 1);

		if (f0 != -1 && f1 != -1)
			edgeArea(i) = (faceArea(f0) + faceArea(f1)) / 2.;
		else if (f0 == -1)
			edgeArea(i) = faceArea(f1) / 2.;
		else
			edgeArea(i) = faceArea(f0) / 2.;
	}
	return edgeArea;
}


Eigen::VectorXd getVertArea(const Eigen::MatrixXd& V, const MeshConnectivity& mesh)
{
	Eigen::VectorXd faceArea = getFaceArea(V, mesh);
	Eigen::VectorXd vertArea;
	vertArea.setZero(V.rows());

	for (int i = 0; i < mesh.nFaces(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int vid = mesh.faceVertex(i, j);
			vertArea(vid) += faceArea(i) / 3.;
		}
	}

	return vertArea;
}
