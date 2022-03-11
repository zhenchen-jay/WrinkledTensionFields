#pragma once

#include <set>
#include <vector>
#include <map>
#include <Eigen/Sparse>

#include "../MeshConnectivity.h"
#include "../MeshGeometry.h"
#include "../Obstacle.h"
#include "../CommonFunctions.h"


class ElasticSetup
{
public:
	ElasticSetup()
	{
		restV.resize(0, 0);
		restF.resize(0, 0);
		clampedDOFs.clear();
		obs.clear();
		restEdgeDOFs.resize(0);

		thickness = 0;
		YoungsModulus = 0;
		PoissonsRatio = 0;
		density = 0;
		penaltyK = 0;
		pressure = 0;

		innerEta = 0;
		restFlat = true;

		maxStepSize = 0;
		perturb = 0;
		gravity.setZero();
		framefreq = 1;
		numInterp = 1;

		abars.clear();
		bbars.clear();
		vertArea.clear();

		isNoeHookean = false;
		tensionField = false;

		bendingType = "";
		sffType = "";
		restMeshPath = "";
		obstaclePath = "";

		initMeshPath = "";
		curMeshPath = "";
		curEdgeDOFsPath = "";

		clampedDOFsPath = "";
		
	}

public:
	Eigen::MatrixXd restV;
	Eigen::MatrixXi restF; // mesh vertices of the original (unstitched) state

	std::vector<Obstacle> obs;

	Eigen::VectorXd restEdgeDOFs;
	Eigen::SparseMatrix<double> laplacian;
	std::map<int, double> clampedDOFs;

	std::string sffType;
	std::shared_ptr<SecondFundamentalFormDiscretization> sff;
	
	double thickness;
	double YoungsModulus;
	double PoissonsRatio;
	double density;
	double penaltyK;
	double innerEta;

	double perturb;

	double pressure;
	bool restFlat;
	bool isNoeHookean;
	bool tensionField;
	std::string bendingType;
	
	double maxStepSize;

	Eigen::Vector3d gravity;
	int framefreq;
	int numInterp; // number of interpolation steps before reaching clamped boundary

	// Derived from the above
	std::vector<Eigen::Matrix2d> abars;
	std::vector<Eigen::Matrix2d> bbars;

	//vert area
	std::vector<double> vertArea;

	std::string restMeshPath, obstaclePath, initMeshPath, curMeshPath, curEdgeDOFsPath, clampedDOFsPath;

public:
	void buildRestFundamentalForms();
	void computeVertArea(const int nverts, const Eigen::MatrixXi& stitchedF);
	void computeLaplacian(const Eigen::MatrixXd restV, const Eigen::MatrixXi restF, const std::vector<Eigen::Vector3i>& bnd_edges, const Eigen::VectorXi& newIndex, const double nverts);
};
