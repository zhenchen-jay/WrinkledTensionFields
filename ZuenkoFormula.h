#pragma once

#include <iostream>
#include "CommonFunctions.h"
#include "MeshGeometry.h"
#include "IntrinsicGeometry.h"
#include "MeshConnectivity.h"


double pointLineSegmentDistance(Eigen::Vector3d p, Eigen::Vector3d bpt, Eigen::Vector3d ept);	// p(2) = bpt(2) = ept(2) = 0

void computeDistanceX(Eigen::MatrixXd baseV, Eigen::MatrixXi baseF, Eigen::MatrixXd restV, Eigen::MatrixXi restF, std::set<int> clampedVerts, Eigen::VectorXd& x);

void computeLambdaAndAmp(Eigen::MatrixXd baseV, Eigen::MatrixXi baseF, Eigen::MatrixXd restV, Eigen::MatrixXi restF, Eigen::VectorXd x, double Ef, double Hf, Eigen::VectorXd& lambda, Eigen::VectorXd& amp);
/* 
	x:		perface value, measure the distance from the centroid to the clamped boundary
	Ef:		Young's modulus
	Hf:		thickness

	lambda = (Ef * Hf / T)^{1/4} * (x / Hf)^{1/2}, T is the tensile strain of the right Cauchy-Green strain, namely, the positive eval of Iu^{-1} (I - Iu)
	amp = lambda / pi (eps - eps^2 / 2)^{1/2}, where eps is the absolute value of Cauchy compression, (L - l) / L. Note that if c to be the negatuve eval of Iu^{-1} (I - Iu), then

	c =  -(l^2 - L^2) / L^2 => eps - eps^2 / 2 = (L^2 - l^2) / (2 L ^2) = - c / 2 => amp = lambda / (sqrt(2) * pi) * sqrt(-c)

*/ 

void getZuenkoPhi(Eigen::MatrixXd baseV, Eigen::MatrixXi baseF, Eigen::MatrixXd restV, Eigen::MatrixXi restF, Eigen::VectorXd lambda, Eigen::VectorXd amp, Eigen::MatrixXd& cutV, Eigen::MatrixXi& cutF, Eigen::VectorXd& phi, Eigen::VectorXd& cutPhi, Eigen::VectorXd& cutAmp, std::set<int> &problemFaces);



void getZuenkoAmpPhi(Eigen::MatrixXd baseV, Eigen::MatrixXi baseF, Eigen::MatrixXd restV, Eigen::MatrixXi restF, std::set<int> clampedVerts, double Ef, double Hf, Eigen::MatrixXd& cutV, Eigen::MatrixXi& cutF, Eigen::VectorXd& phi, Eigen::VectorXd& amp, Eigen::VectorXd& dphi, Eigen::VectorXd& cutPhi, Eigen::VectorXd& cutAmp, std::set<int>& problemFaces);
/* compute amp and phi using the formula proposed in Zuenko's paper

Inputs:
	(baseV, baseF):		the base mesh (vertices, faces)
	(restV, restF):		the rest mesh (vertices, faces)
	clampedVerts:		the clampled vertices id, note that all the id are based on base mesh
	Ef:					Young's modulus
	Hf:					thickness

Outputs:
	(cutV, cutF):		the cutted base mesh (vertices, faces), due to the singularities or pure tension faces of the vector field
	(amp, phi, dphi):	the recovered phi, amp and dphi on the uncutted base mesh, note that there will be some interger of 2pi jump in the phi along the cuts (per vertex value)
	(cutAmp, cutPhi):	the recovered phi and amp on the cutted base mesh (per vertex value)
	problemFaces:		a set containing all the problem faces, namely, the face amplitude is 0.
	

*/