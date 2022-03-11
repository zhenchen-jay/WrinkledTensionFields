#include <igl/boundary_loop.h>
#include <igl/boundary_facets.h>
#include <igl/on_boundary.h>

#include "ZuenkoFormula.h"
#include "WTFShell/PhiEstimate.h"

static const double PI = 3.1415926535898;

double pointLineSegmentDistance(Eigen::Vector3d p, Eigen::Vector3d bpt, Eigen::Vector3d ept)
{
	Eigen::Vector3d pe, be, pb;
	pe = ept - p;
	pb = bpt - p;
	be = ept - bpt;

	double c = pe.dot(be) / be.dot(be);

	if (c >= 0 && c <= 1)
	{
		Eigen::Vector3d q = c * bpt + (1 - c) * ept;
		return (p - q).norm();
	}

	else
		return std::min((p - bpt).norm(), (p - ept).norm());
}

void computeDistanceX(Eigen::MatrixXd baseV, Eigen::MatrixXi baseF, Eigen::MatrixXd restV, Eigen::MatrixXi restF, std::set<int> clampedVerts, Eigen::VectorXd& x)
{
	std::cout << "clamped vertices size: " <<clampedVerts.size() << std::endl;

	MeshConnectivity restMesh(restF);
	MeshConnectivity baseMesh(baseF);

	MeshGeometry baseGeo(baseV, baseMesh);
	MeshGeometry restGeo(restV, restMesh);

	int nfaces = baseF.rows();

	x.resize(nfaces);
	x.setZero();

	Eigen::MatrixXd centroid(nfaces, 3);
	centroid.setZero();

	for (int i = 0; i < nfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			centroid.row(i) += restV.row(restMesh.faceVertex(i, j)) / 3.0;
		}
	}

	std::vector<int> clampedBndEdges;
	for (int i = 0; i < restMesh.nEdges(); i++)
	{
		if (restMesh.edgeFace(i, 0) != -1 && restMesh.edgeFace(i, 1) != -1) // not on the boundary
			continue;
		// since the clamped vertex id is the baseV vid, we need to convert to get the correct vid
		int fid = restMesh.edgeFace(i, 0) != -1 ? restMesh.edgeFace(i, 0) : restMesh.edgeFace(i, 1);

		int v0 = restMesh.edgeVertex(i, 0);
		int v1 = restMesh.edgeVertex(i, 1);

		int vid0 = -1;
		int vid1 = -1;
		for (int j = 0; j < 3; j++)
		{
			if (restMesh.faceVertex(fid, j) == v0)
				vid0 = j;
			if (restMesh.faceVertex(fid, j) == v1)
				vid1 = j;
		}
		int baseV0 = baseMesh.faceVertex(fid, vid0);
		int baseV1 = baseMesh.faceVertex(fid, vid1);

		if (clampedVerts.count(baseV0) && clampedVerts.count(baseV1)) // this edge is clamped
		{
			clampedBndEdges.push_back(i);
		}
	}

	for (int i = 0; i < nfaces; i++)
	{
		double minDis = std::numeric_limits<double>::infinity();
		for (int j = 0; j < clampedBndEdges.size(); j++)
		{
			int eid = clampedBndEdges[j];
			int v0 = restMesh.edgeVertex(eid, 0);
			int v1 = restMesh.edgeVertex(eid, 1);
			double d = pointLineSegmentDistance(centroid.row(i), restV.row(v0), restV.row(v1));

			minDis = std::min(minDis, d);
		}
		x(i) = minDis;
	}
	std::cout << "distance  (max, min): " << x.maxCoeff() << ", " << x.minCoeff() << std::endl;
}

void computeLambdaAndAmp(Eigen::MatrixXd baseV, Eigen::MatrixXi baseF, Eigen::MatrixXd restV, Eigen::MatrixXi restF, Eigen::VectorXd x, double Ef, double Hf, Eigen::VectorXd& lambda, Eigen::VectorXd& amp)
{
	MeshConnectivity restMesh(restF);
	MeshConnectivity baseMesh(baseF);

	MeshGeometry baseGeo(baseV, baseMesh);
	MeshGeometry restGeo(restV, restMesh);

	int nfaces = baseF.rows();

	lambda.resize(nfaces);
	lambda.setZero();

	Eigen::VectorXd perFaceAmp(nfaces);
	perFaceAmp.setZero();


	for (int i = 0; i < nfaces; i++)
	{
		Eigen::Matrix2d abar = restGeo.Bs[i].transpose() * restGeo.Bs[i];
		Eigen::Matrix2d a = baseGeo.Bs[i].transpose() * baseGeo.Bs[i];
		Eigen::Matrix2d strain = abar.inverse() * (a - abar);
		Eigen::EigenSolver<Eigen::Matrix2d> solver(strain);
		
		Eigen::Vector2d evals = solver.eigenvalues().real();
		Eigen::Matrix2d evecs = solver.eigenvectors().real();


		double compression = evals(0);
		double tensile = evals(1);

		if (compression > tensile)
			std::swap(compression, tensile);

		if (compression >= 0 || tensile <= 0)
			continue;

		lambda(i) = Hf * std::pow(Ef * Hf / (1 + tensile), 1.0 / 4.0) * std::sqrt(x(i) / Hf);
		
		perFaceAmp(i) = lambda(i) / PI * sqrt(-compression / 2);
	}

	amp.resize(baseV.rows());
	amp.setZero();

	Eigen::VectorXd vertDegree(baseV.rows());
	vertDegree.setZero();

	for (int i = 0; i < nfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int vid = baseMesh.faceVertex(i, j);
			vertDegree(vid) += 1.0;
			amp(vid) += perFaceAmp(i);
		}
	}

	for (int i = 0; i < amp.size(); i++)
	{
		amp(i) /= vertDegree(i);
	}

	std::cout <<"E: " << Ef << ", Hf: " << Hf <<", lambda (min, max): " << lambda.minCoeff() << ", " << lambda.maxCoeff() << ", amp (min, max): " << amp.minCoeff() << ", " << amp.maxCoeff() << std::endl;
}

void getZuenkoPhi(Eigen::MatrixXd baseV, Eigen::MatrixXi baseF, Eigen::MatrixXd restV, Eigen::MatrixXi restF, Eigen::VectorXd lambda, Eigen::VectorXd amp, Eigen::MatrixXd& cutV, Eigen::MatrixXi& cutF, Eigen::VectorXd& phi, Eigen::MatrixXd& w, Eigen::VectorXd& dphi, Eigen::VectorXd& cutPhi, Eigen::VectorXd &cutAmp, std::set<int>& problemFaces)
{
	std::vector<Eigen::Matrix2d> abars;

	MeshConnectivity restMesh(restF);
	MeshConnectivity baseMesh(baseF);

	MeshGeometry baseGeo(baseV, baseMesh);
	MeshGeometry restGeo(restV, restMesh);

	int nfaces = baseF.rows();

	Eigen::MatrixXd wguess(nfaces, 2);
	wguess.setZero();

	
	for (int i = 0; i < nfaces; i++)
	{
		Eigen::Matrix2d abar = restGeo.Bs[i].transpose() * restGeo.Bs[i];
		abars.push_back(abar);
	}

	locatePotentialPureTensionFaces(abars, baseV, baseF, problemFaces);


	for (int i = 0; i < nfaces; i++)
	{
		Eigen::Matrix2d abar = abars[i];
		Eigen::Matrix2d a = baseGeo.Bs[i].transpose() * baseGeo.Bs[i];
		Eigen::Matrix2d diff = a - abar;
		Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix2d> solver(diff, abar);
		Eigen::Vector2d evec = solver.eigenvectors().col(0);

		double compression = solver.eigenvalues()(0);
		double tensile = solver.eigenvalues()(1);

		if (problemFaces.find(i) != problemFaces.end() || lambda(i) == 0)
			continue;

		wguess.row(i) = 2 * PI / lambda(i) * (abar * evec).transpose();	// convert to one-form
		
	}
	std::cout << "||w|| = " << wguess.norm() << std::endl;
	combField(baseF, abars, NULL, wguess, w);

	roundPhiFromOmegaCutbyTension(baseV, baseF, cutV, cutF, abars, amp, w, phi, cutPhi, cutAmp, problemFaces);

	MeshConnectivity cutMesh(cutF);
	Eigen::VectorXd cutDphi(cutMesh.nEdges());
	for (int i = 0; i < cutMesh.nEdges(); i++)
	{
		int v0 = cutMesh.edgeVertex(i, 0);
		int v1 = cutMesh.edgeVertex(i, 1);

		if (v0 < v1)
			cutDphi(i) = cutPhi(v1) - cutPhi(v0);
		else
			cutDphi(i) = cutPhi(v0) - cutPhi(v1);


	}

	int nedges = baseMesh.nEdges();
	dphi.resize(nedges);
	dphi.setZero();

	Eigen::VectorXi isVisited(baseMesh.nEdges());
	isVisited.setZero();

	for (int f = 0; f < nfaces; f++)
	{
		if (problemFaces.find(f) != problemFaces.end())
			continue;
		for (int e = 0; e < 3; e++)
		{
			int eid = baseMesh.faceEdge(f, e);
			if (isVisited(eid))
				continue;
			isVisited(eid) = 1;

			int vid0 = baseMesh.edgeVertex(eid, 0);
			int vid1 = baseMesh.edgeVertex(eid, 1);
			Eigen::Vector2d phiVals;
			phiVals << std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();

			for (int v = 0; v < 3; v++)
			{
				int cutVid = cutMesh.faceVertex(f, v);

				if ((cutV.row(cutVid) - baseV.row(vid0)).norm() < 1e-6)
					phiVals(0) = cutPhi(cutVid);

				if ((cutV.row(cutVid) - baseV.row(vid1)).norm() < 1e-6)
					phiVals(1) = cutPhi(cutVid);

			}

			if (phiVals(0) == std::numeric_limits<double>::infinity() || phiVals(1) == std::numeric_limits<double>::infinity())
				std::cout << "error!" << std::endl;


			if (vid0 < vid1)
				dphi(eid) = phiVals(1) - phiVals(0);
			else
				dphi(eid) = phiVals(1) - phiVals(0);

		}
	}
	//faceDPhi2EdgeDPhi(w, baseF, dphi);
}

void getZuenkoAmpPhi(Eigen::MatrixXd baseV, Eigen::MatrixXi baseF, Eigen::MatrixXd restV, Eigen::MatrixXi restF, std::set<int> clampedVerts, double Ef, double Hf, Eigen::MatrixXd& cutV, Eigen::MatrixXi& cutF, Eigen::VectorXd& phi, Eigen::VectorXd& amp, Eigen::VectorXd& dphi, Eigen::VectorXd& cutPhi, Eigen::VectorXd& cutAmp, std::set<int>& problemFaces)
{
	Eigen::VectorXd dist;
	computeDistanceX(baseV, baseF, restV, restF, clampedVerts, dist);

	Eigen::VectorXd lambda;
	computeLambdaAndAmp(baseV, baseF, restV, restF, dist, Ef, Hf, lambda, amp);

	Eigen::MatrixXd w;

	getZuenkoPhi(baseV, baseF, restV, restF, lambda, amp, cutV, cutF, phi, w, dphi, cutPhi, cutAmp, problemFaces);
	
	for(auto & v : clampedVerts)
		amp(v) = 0;
}