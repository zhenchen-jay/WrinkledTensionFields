#include"ElasticState.h"

ElasticState::ElasticState(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const SecondFundamentalFormDiscretization& sff)
{
	initialGuess = V;
	mesh = MeshConnectivity(F);
	curPos = V;

	sff.initializeExtraDOFs(initialEdgeDOFs, mesh, V);
	curEdgeDOFs = initialEdgeDOFs;
}