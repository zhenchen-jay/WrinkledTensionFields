#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../CommonFunctions.h"
#include "../Obstacle.h"

namespace WTF
{
	class WTFSetup
	{
	public:
		WTFSetup()
		{
			restV.resize(0, 0);
			restF.resize(0, 0);
			clampedDOFs.clear();
			obs.clear();
			restEdgeDOFs.resize(0);

			thickness = 0;
			YoungsModulus = 0;
			PoissonsRatio = 0;
			nasoqEps = 0;

			quadPoints.clear();
			abars.clear();
			bbars.clear(); 
			
			clampedChosenVerts = true;
			restFlat = true;
			quadNum = 0;
			sffType = "";
			restMeshPath = "";
			obstaclePath = "";
			baseMeshPath = "";
			clampedDOFsPath = "";
			phiPath = "";
			dphiPath = "";
			ampPath = "";

		}

	public:
		// Core data structures
		Eigen::MatrixXd restV;
		Eigen::MatrixXi restF; // mesh vertices of the original (unstitched) state
		std::vector<Obstacle> obs;

		std::map<int, double> clampedDOFs;
		double thickness;
		double YoungsModulus;
		double PoissonsRatio;
		double nasoqEps;

		int quadNum;
		std::string sffType;
		
		bool restFlat;
		bool clampedChosenVerts;
		std::shared_ptr<SecondFundamentalFormDiscretization> sff;
		Eigen::VectorXd restEdgeDOFs;

		std::vector<QuadraturePoints> quadPoints;
		
		
		// Derived from the above
		std::vector<Eigen::Matrix2d> abars;
		std::vector<Eigen::Matrix2d> bbars;

	public:	// file pathes
		std::string restMeshPath, obstaclePath, baseMeshPath, clampedDOFsPath;
		std::string phiPath, dphiPath, ampPath;

		void buildQuadraturePoints();
		void buildRestFundamentalForms();
	};
}