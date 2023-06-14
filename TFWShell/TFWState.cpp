#include "TFWState.h"
#include "PhiEstimate.h"
#include "MeshUpsampling.h"

#include "../CommonFunctions.h"

using namespace TFW;
void TFWState::getWrinkleMesh(const TFWSetup& setup, int upsamplingTimes, bool isUseV1Term, bool isUseV2Term, RoundingType roundingType, bool isFixBnd)
{
	Eigen::MatrixXd cutV;
	Eigen::MatrixXi cutF;
	Eigen::VectorXd cutPhi, cutDphi, cutAmp;
	std::set<int> clampedVerts;
	clampedVerts.clear();

	if (setup.clampedChosenVerts)
	{
		for (auto& it : setup.clampedDOFs)
		{
			int vid = it.first / 3;
			if (clampedVerts.count(vid) == 0)
				clampedVerts.insert(vid);
		}
	}

	if(roundingType == CWFRound)
	{
		wrinkledMeshUpsampling(basePos, baseMesh.faces(), setup.restV, setup.restF, amplitude, dphi, phi, &wrinkledPos, &wrinkledF, nullptr, nullptr, nullptr, nullptr, *(setup.sff), setup.YoungsModulus, setup.PoissonsRatio, upsamplingTimes, isFixBnd, isUseV1Term, isUseV2Term);
	}
	else
	{
		roundPhiFromDphiCutbyTension(basePos, baseMesh.faces(), cutV, cutF, setup.abars, amplitude, dphi, roundingType, phi, cutPhi, cutAmp, tensionFaces);
		wrinkledMeshUpsamplingUncut(basePos, baseMesh.faces(), setup.restV, setup.restF, cutV, cutF, tensionFaces, clampedVerts,
		                            &wrinkledPos, &wrinkledF, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
		                            cutAmp, cutPhi, *(setup.sff), setup.YoungsModulus, setup.PoissonsRatio, upsamplingTimes, Loop, isUseV1Term, isUseV2Term);
	}



}


void TFWState::reinitializeWrinkleVaribles(const TFWSetup& setup, RoundingType roundingType)
{
	std::set<int> clampedVerts;
	clampedVerts.clear();
	std::map<int, double>::const_iterator it;
	for (it = setup.clampedDOFs.begin(); it != setup.clampedDOFs.end(); it++)
	{
		int vid = it->first / 3;
		if (clampedVerts.find(vid) == clampedVerts.end() && setup.clampedChosenVerts)
			clampedVerts.insert(vid);
	}

	std::cout << "Reinitialize amp and dphi." << std::endl;
	estimateWrinkleVariablesFromStrainCutbyTension(setup.abars, basePos, baseMesh.faces(), clampedVerts, 0.01 * (basePos.maxCoeff() - basePos.minCoeff()), amplitude, phi, dphi, tensionFaces, roundingType);
	dualAmp.resize(0);
	dualDphi.resize(0);
}
