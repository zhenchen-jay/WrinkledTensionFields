#ifndef COMMONFUNCTIONS_H
#define COMMONFUNCTIONS_H

#include <set>

#include <iostream>
#include <fstream>

#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include "SecondFundamentalForm/SecondFundamentalFormDiscretization.h"
#include "SecondFundamentalForm/MidedgeAngleSinFormulation.h"
#include "SecondFundamentalForm/MidedgeAngleTanFormulation.h"
#include "SecondFundamentalForm/MidedgeAverageFormulation.h"

struct QuadraturePoints
{
	double u;
	double v;
	double weight;
};


enum StretchingType
{
	StVK = 0,
	tensionField = 1,
	NeoHookean = 2
};

enum BendingType
{
	elasticBending = 0,
	quadraticBending = 1,
	noBending = 2
};

enum SolverType
{
	Lbfgs = 0,
	Newton = 1,
	ActiveSet = 2
};

enum SFFType
{
	MidedgeAverage = 0,
	MidedgeSin = 1,
	MidedgeTan = 2
};

enum SubdivisionType
{
	Midpoint = 0,
	Loop = 1
};

enum RoundingType
{
	ComisoRound = 0,
	CWFRound = 1
};

struct TimeCost
{
	double gradTime = 0;
	double hessTime = 0;
	double solverTime = 0;
	double collisionDectionTime = 0;
	double lineSearchTime = 0;
	double updateTime = 0;
	double convergenceCheckTime = 0;
	double savingTime = 0;
	double unconstrainedLLTTime = 0;

	double totalTime() 
	{
		return gradTime + hessTime + solverTime + collisionDectionTime + lineSearchTime + updateTime + convergenceCheckTime; // exclude the saving time
	}
};

Eigen::MatrixXd lowRankApprox(Eigen::MatrixXd A);      // semi-positive projection of a symmetric matrix A


double cotan(const Eigen::Vector3d v0, const Eigen::Vector3d v1, const Eigen::Vector3d v2);

void locatePotentialPureTensionFaces(const std::vector<Eigen::Matrix2d>& abars, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, std::set<int>& potentialPureTensionFaces);
void getPureTensionVertsEdges(const std::set<int>& potentialPureTensionFaces, const Eigen::MatrixXi& F, std::set<int>* pureTensionEdges, std::set<int> *pureTensionVerts);

Eigen::MatrixXd nullspaceExtraction(const Eigen::SparseMatrix<double> A);

void trivialOffset(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd &offsettedV, double d); // offfset V along its normal direction by a given distance

void rigidBodyAlignment(const Eigen::MatrixXd& tarPos, const MeshConnectivity& mesh, const Eigen::MatrixXd &pos, Eigen::Matrix3d &R, Eigen::Vector3d &t); // implement a really naive iterative closest point (ICP) algorithm for rigid body alignment.

void matToVec(const Eigen::MatrixXd& mat, Eigen::VectorXd& vec);
void vecToMat(const Eigen::VectorXd& vec, Eigen::MatrixXd& mat);

bool mkdir(const std::string& foldername);

Eigen::MatrixXd intrinsicHalfEdgeVec2VertexVec(const Eigen::MatrixXd& v, const Eigen::MatrixXd& pos, const MeshConnectivity& mesh);
Eigen::MatrixXd intrinsicHalfEdgeVec2FaceVec(const Eigen::MatrixXd& w, const Eigen::MatrixXd& pos, const MeshConnectivity& mesh);
Eigen::MatrixXd intrinsicEdgeVec2FaceVec(const Eigen::VectorXd& w, const Eigen::MatrixXd& pos, const MeshConnectivity& mesh);

Eigen::VectorXd faceVec2IntrinsicEdgeVec(const Eigen::MatrixXd& v, const Eigen::MatrixXd& pos, const MeshConnectivity& mesh);
Eigen::MatrixXd intrinsicEdgeVec2FaceVec(const Eigen::VectorXd& v, const Eigen::MatrixXd& pos, const MeshConnectivity& mesh);
Eigen::MatrixXd intrinsicHalfEdgeVec2FaceVec(const Eigen::MatrixXd& v, const Eigen::MatrixXd& pos, const MeshConnectivity& mesh);

void computeBaryGradient(const Eigen::Vector3d& P0, const Eigen::Vector3d& P1, const Eigen::Vector3d& P2, const Eigen::Vector3d& bary, Eigen::Matrix3d& baryGrad);

Eigen::VectorXd getFaceArea(const Eigen::MatrixXd& V, const MeshConnectivity& mesh);
Eigen::VectorXd getEdgeArea(const Eigen::MatrixXd& V, const MeshConnectivity& mesh);
Eigen::VectorXd getVertArea(const Eigen::MatrixXd& V, const MeshConnectivity& mesh);

#endif