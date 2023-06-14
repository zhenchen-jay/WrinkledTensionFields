#include "ComplexLoopZuenko.h"
#include "../CommonFunctions.h"
#include <iostream>
#include <cassert>
#include <memory>


std::complex<double> ComplexLoopZuenko::interpZ(const std::vector<std::complex<double>>& zList, const std::vector<Eigen::Vector3d>& gradThetaList, std::vector<double>& coords, const std::vector<Eigen::Vector3d>& pList)
{
	Eigen::Vector3d P = Eigen::Vector3d::Zero();
	
	int n = pList.size();
	for (int i = 0; i < n; i++)
	{
		P += coords[i] * pList[i];
	}
	Eigen::Vector3d gradThetaAve = Eigen::Vector3d::Zero();
	for (int i = 0; i < n; i++)
	{
		gradThetaAve += coords[i] * gradThetaList[i];
	}

	std::complex<double> z = 0;

	std::complex<double> I = std::complex<double>(0, 1);

	std::vector<double> deltaThetaList(n);
	for (int i = 0; i < n; i++)
	{
		//deltaThetaList[i] = (P - pList[i]).dot(gradThetaList[i]);
		deltaThetaList[i] = (P - pList[i]).dot(gradThetaAve);
	}

	Eigen::VectorXcd dzdalpha;
	dzdalpha.setZero(n);

	for (int i = 0; i < n; i++)
	{
		z += coords[i] * zList[i] * std::complex<double>(std::cos(deltaThetaList[i]), std::sin(deltaThetaList[i]));
	}

	return z;
}

std::complex<double> ComplexLoopZuenko::computeZandGradZ(const Eigen::VectorXd& omega, const std::vector<std::complex<double>>& zvals, int fid, const Eigen::Vector3d& bary, Eigen::Vector3cd* gradz)
{
	Eigen::Matrix3d gradBary;
	Eigen::Vector3d P0 = _mesh.GetVertPos(_mesh.GetFaceVerts(fid)[0]);
	Eigen::Vector3d P1 = _mesh.GetVertPos(_mesh.GetFaceVerts(fid)[1]);
	Eigen::Vector3d P2 = _mesh.GetVertPos(_mesh.GetFaceVerts(fid)[2]);
	
	computeBaryGradient(P0, P1, P2, bary, gradBary);

	std::complex<double> z = 0;

	if (gradz)
	{
		gradz->setZero();
	}
	
	
	std::complex<double> tau = std::complex<double>(0, 1);

	for (int i = 0; i < 3; i++)
	{
		double alphai = bary(i);
		int vid = _mesh.GetFaceVerts(fid)[i];
		int eid0 = _mesh.GetFaceEdges(fid)[i];
		int eid1 = _mesh.GetFaceEdges(fid)[(i + 2) % 3];

		double w1 = omega(eid0);
		double w2 = omega(eid1);

		if (vid == _mesh.GetEdgeVerts(eid0)[1])
			w1 *= -1;
		if (vid == _mesh.GetEdgeVerts(eid1)[1])
			w2 *= -1;

		double deltaTheta = w1 * bary((i + 1) % 3) + w2 * bary((i + 2) % 3);
		std::complex<double> expip = std::complex<double>(std::cos(deltaTheta), std::sin(deltaTheta));

		z += alphai * zvals[vid] * expip;

		if(gradz)
			(*gradz) += zvals[vid] * expip * (gradBary.row(i) + tau * alphai * w1 * gradBary.row((i + 1) % 3) + tau * alphai * w2 * gradBary.row((i + 2) % 3));
	}
	return z;
}

Eigen::Vector3d ComplexLoopZuenko::computeGradThetaFromOmegaPerface(const Eigen::VectorXd& omega, int fid, int vInF)
{
    int vid = _mesh.GetFaceVerts(fid)[vInF];
    int eid0 = _mesh.GetFaceEdges(fid)[vInF];
    int eid1 = _mesh.GetFaceEdges(fid)[(vInF + 2) % 3];
    Eigen::Vector3d r0 = _mesh.GetVertPos(_mesh.GetEdgeVerts(eid0)[1]) - _mesh.GetVertPos(_mesh.GetEdgeVerts(eid0)[0]);
    Eigen::Vector3d r1 = _mesh.GetVertPos(_mesh.GetEdgeVerts(eid1)[1]) - _mesh.GetVertPos(_mesh.GetEdgeVerts(eid1)[0]);

    Eigen::Matrix2d Iinv, I;
    I << r0.dot(r0), r0.dot(r1), r1.dot(r0), r1.dot(r1);
    Iinv = I.inverse();

    Eigen::Vector2d rhs;
    double w1 = omega(eid0);
    double w2 = omega(eid1);
    rhs << w1, w2;

    Eigen::Vector2d u = Iinv * rhs;
    return u[0] * r0 + u[1] * r1;
}

Eigen::Vector3d ComplexLoopZuenko::computeBaryGradThetaFromOmegaPerface(const Eigen::VectorXd& omega, int fid, const Eigen::Vector3d& bary)
{
    Eigen::Vector3d gradTheta = Eigen::Vector3d::Zero();
    for(int i = 0; i < 3; i++)
    {
        gradTheta += bary[i] * computeGradThetaFromOmegaPerface(omega, fid, i);
    }
    return gradTheta;
}

void ComplexLoopZuenko::updateLoopedZvals(const Eigen::VectorXd& omega, const std::vector<std::complex<double>>& zvals, std::vector<std::complex<double>>& upZvals)
{
	int V = _mesh.GetVertCount();
	int E = _mesh.GetEdgeCount();
	
	upZvals.resize(V + E);
	// Even verts
	for (int vi = 0; vi < V; ++vi)
	{
		if (_mesh.IsVertBoundary(vi))
		{
			if(_isFixBnd)
				upZvals[vi] = zvals[vi];
			else
			{
				std::vector<int> boundary(2);
				boundary[0] = _mesh.GetVertEdges(vi).front();
				boundary[1] = _mesh.GetVertEdges(vi).back();

				std::vector<std::complex<double>> zp(2);
				std::vector<Eigen::Vector3d> gradthetap(2);
				std::vector<double> coords = { 1. / 2, 1. / 2 };
				std::vector<Eigen::Vector3d> pList(2);

				for (int j = 0; j < boundary.size(); ++j)
				{
					int edge = boundary[j];
					int face = _mesh.GetEdgeFaces(edge)[0];
					int viInface = _mesh.GetVertIndexInFace(face, vi);

					int viInEdge = _mesh.GetVertIndexInEdge(edge, vi);
					int vj = _mesh.GetEdgeVerts(edge)[(viInEdge + 1) % 2];

					int vjInface = _mesh.GetVertIndexInFace(face, vj);

					Eigen::Vector3d bary = Eigen::Vector3d::Zero();
					bary(viInface) = 3. / 4;
					bary(vjInface) = 1. / 4;

                    pList[j] = 3. / 4 * _mesh.GetVertPos(vi) + 1. / 4 * _mesh.GetVertPos(vj);
                    // grad from vi
                    gradthetap[j] = computeBaryGradThetaFromOmegaPerface(omega, face, bary);
                    zp[j] = computeZandGradZ(omega, zvals, face, bary, NULL);


                    // this is really unstable
//					Eigen::Vector3cd gradZ;
//					zp[j] = computeZandGradZ(omega, zvals, face, bary, &(gradZ));
//
//
//					gradthetap[j] = (std::conj(zp[j]) * gradZ).imag();
//
//					if (std::abs(zp[j]))
//						gradthetap[j] = gradthetap[j] / (std::abs(zp[j]) * std::abs(zp[j]));

				}
				upZvals[vi] = interpZ(zp, gradthetap, coords, pList);
			}
			
			
		}
		else
		{
			const std::vector<int>& vFaces = _mesh.GetVertFaces(vi);
			int nNeiFaces = vFaces.size();

			// Fig5 left [Wang et al. 2006]
			Scalar alpha = 0.375;
			if (nNeiFaces == 3) alpha /= 2;
			else                    alpha /= nNeiFaces;

			double beta = nNeiFaces / 2. * alpha;

			std::vector<std::complex<double>> zp(nNeiFaces);
			std::vector<Eigen::Vector3d> gradthetap(nNeiFaces);
			std::vector<double> coords;
			coords.resize(nNeiFaces, 1. / nNeiFaces);
			std::vector<Eigen::Vector3d> pList(nNeiFaces);

			for (int k = 0; k < nNeiFaces; ++k)
			{
				int face = vFaces[k];
				int viInface = _mesh.GetVertIndexInFace(face, vi);
				Eigen::Vector3d bary;
				bary.setConstant(beta);
				bary(viInface) = 1 - 2 * beta;

				pList[k] = Eigen::Vector3d::Zero();
				for (int i = 0; i < 3; i++)
				{
					pList[k] += bary(i) * _mesh.GetVertPos(_mesh.GetFaceVerts(face)[i]);
				}

                zp[k] = computeZandGradZ(omega, zvals, face, bary, NULL);
                gradthetap[k] = computeBaryGradThetaFromOmegaPerface(omega, face, bary);
            

//				Eigen::Vector3cd gradZ;
//				zp[k] = computeZandGradZ(omega, zvals, face, bary, &(gradZ));
//
//				gradthetap[k] = (std::conj(zp[k]) * gradZ).imag();
//
//				if(std::abs(zp[k]))
//					gradthetap[k] = gradthetap[k] / (std::abs(zp[k]) * std::abs(zp[k]));
			}
			upZvals[vi] = interpZ(zp, gradthetap, coords, pList);
		}
	}

	// Odd verts
	for (int edge = 0; edge < E; ++edge)
	{
		int row = edge + V;
		if (_mesh.IsEdgeBoundary(edge))
		{
			int face = _mesh.GetEdgeFaces(edge)[0];
			int eindexFace = _mesh.GetEdgeIndexInFace(face, edge);
			Eigen::Vector3d bary;
			bary.setConstant(0.5);
			bary((eindexFace + 2) % 3) = 0;
            
			upZvals[row] = computeZandGradZ(omega, zvals, face, bary, NULL);
		}
		else
		{
			
			std::vector<std::complex<double>> zp(2);
			std::vector<Eigen::Vector3d> gradthetap(2);
			std::vector<double> coords = { 1. / 2, 1. / 2 };
			std::vector<Eigen::Vector3d> pList(2);

			for (int j = 0; j < 2; ++j)
			{
				int face = _mesh.GetEdgeFaces(edge)[j];
				int offset = _mesh.GetEdgeIndexInFace(face, edge);

				Eigen::Vector3d bary;
				bary.setConstant(3. / 8.);
				bary((offset + 2) % 3) = 0.25;

				pList[j] = Eigen::Vector3d::Zero();
				for (int i = 0; i < 3; i++)
				{
					pList[j] += bary(i) * _mesh.GetVertPos(_mesh.GetFaceVerts(face)[i]);
				}

                zp[j] = computeZandGradZ(omega, zvals, face, bary, NULL);
                gradthetap[j] = computeBaryGradThetaFromOmegaPerface(omega, face, bary);
             
//				Eigen::Vector3cd gradZ;
//				zp[j] = computeZandGradZ(omega, zvals, face, bary, &(gradZ));
//
//				gradthetap[j] = (std::conj(zp[j]) * gradZ).imag();
//
//				if(std::abs(zp[j]))
//					gradthetap[j] = gradthetap[j] / (std::abs(zp[j]) * std::abs(zp[j]));

			}

			upZvals[row] = interpZ(zp, gradthetap, coords, pList);
		}
	}
}

void ComplexLoopZuenko::Subdivide(const Eigen::VectorXd& omega, const std::vector<std::complex<double>>& zvals, Eigen::VectorXd& omegaNew, std::vector<std::complex<double>>& upZvals, int level)
{
	
	int nverts = _mesh.GetVertCount();
	omegaNew = omega;
	upZvals = zvals;


	MatrixX X;
	_mesh.GetPos(X);

	Eigen::VectorXd amp(nverts);
	
	for (int i = 0; i < nverts; i++)
	{
		amp(i) = std::abs(zvals[i]);
	}

	for (int l = 0; l < level; ++l) 
	{
		SparseMatrixX tmpS0, tmpS1;
		BuildS0(tmpS0);
		BuildS1(tmpS1);

		X = tmpS0 * X;
		amp = tmpS0 * amp;

		std::vector<std::complex<double>> upZvalsNew;
		//std::vector<Eigen::Vector3cd> upGradZvals;

		updateLoopedZvals(omegaNew, upZvals, upZvalsNew);

		//updateLoopedZvalsNew(meshNew, omegaNew, upZvals, upZvalsNew);
		omegaNew = tmpS1 * omegaNew;

		std::vector<Vector3> points;
		ConvertToVector3(X, points);

		std::vector< std::vector<int> > edgeToVert;
		GetSubdividedEdges(edgeToVert);

		std::vector< std::vector<int> > faceToVert;
		GetSubdividedFaces(faceToVert);

		_mesh.Populate(points, faceToVert, edgeToVert);

		upZvals.swap(upZvalsNew);
	}
}
