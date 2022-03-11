#pragma once

#include <Eigen/Core>
#include <Eigen/SPQRSupport>
#include <set>
#include <igl/principal_curvature.h>

#include "../MeshConnectivity.h"
#include "../CommonFunctions.h"

struct ElasticState
{
	ElasticState()
	{
		initialGuess.resize(0, 0);
		initialEdgeDOFs.resize(0);

		curPos = initialGuess;
		curEdgeDOFs = initialEdgeDOFs;
	}

	ElasticState(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const SecondFundamentalFormDiscretization& sff);

	MeshConnectivity mesh;

	Eigen::MatrixXd initialGuess;
	Eigen::VectorXd initialEdgeDOFs;

	Eigen::MatrixXd curPos;
	Eigen::VectorXd curEdgeDOFs;
	
	std::set<int> collidingPoints;
	std::vector<double> mass;

	void resetState()
	{
		curPos = initialGuess;
		curEdgeDOFs = initialEdgeDOFs;
	}
};
