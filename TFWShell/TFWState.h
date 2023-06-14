#pragma once

#include <Eigen/Core>
#include <Eigen/SPQRSupport>
#include <set>
#include <igl/principal_curvature.h>

#include "../MeshLib/MeshConnectivity.h"

#include "TFWSetup.h"

namespace TFW
{
	class TFWState
	{
	public:
		TFWState()
		{
			basePos.resize(0, 0);
			baseEdgeDOFs.resize(0);
			phi.resize(0);
			dphi.resize(0);
			amplitude.resize(0);
			dualAmp.resize(0);
			dualDphi.resize(0);

			tensionFaces.clear();
			PD1.resize(0, 3);
			PD2.resize(0, 3);
			PV1.resize(0);
			PV2.resize(0);

			wrinkledPos = basePos;
			wrinkledF.resize(0, 3);

		}

		Eigen::MatrixXd basePos;
		MeshConnectivity baseMesh;
		Eigen::VectorXd baseEdgeDOFs;
		Eigen::VectorXd phi;
		Eigen::VectorXd dphi;
		Eigen::VectorXd amplitude;

		Eigen::VectorXd dualDphi;
		Eigen::VectorXd dualAmp;

		std::set<int> tensionFaces;

		//Base mesh curvature
		Eigen::MatrixXd PD1;
		Eigen::MatrixXd PD2;
		Eigen::VectorXd PV1;
		Eigen::VectorXd PV2;

		// wrinkle mesh
		Eigen::MatrixXd wrinkledPos;
		Eigen::MatrixXi wrinkledF;

		void computeBaseCurvature()
		{
			if (baseMesh.nFaces() > 1)
				igl::principal_curvature(basePos, baseMesh.faces(), PD1, PD2, PV1, PV2);
		}

		void getWrinkleMesh(const TFWSetup& setup, int upsamplingTimes = 2, bool isUseV1Term = false, bool isUseV2Term = true, RoundingType roundingType = CWFRound, bool isFixBnd = false);
		void reinitializeWrinkleVaribles(const TFWSetup& setup, RoundingType roundingType = CWFRound);

	};
}
