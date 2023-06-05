#include <igl/writeOBJ.h>
#include <igl/null.h>
#include <memory>
#include <set>
#include "TFWShell.h"
#include "TFWModel.h"
#include "MeshUpsampling.h"
#include "PhiEstimate.h"

void TFWModel::initialization(const TFWSetup setup, const TFWState state, std::map<int, double>* clampedDOFs, std::string filePath, bool isUsePosHess, bool isParallel)
{
	_setup = setup;
	//_prevPos = setup.initialPos;
	_state = state;
	setProjM(clampedDOFs);
	_isUsePosHess = isUsePosHess;
	_isParallel = isParallel;
	_filePrefix = filePath;

	if (_isParallel)
		std::cout << "Use TBB for parallel computing the energy per face" << std::endl;
	else
		std::cout << "Sequential computing the energy per face" << std::endl;
}


void TFWModel::setProjM(std::map<int, double>* clampedDOFs)
{
	_freeAmp = 0;
	_freePhi = 0;
	_freeDphi = 0;
	MeshConnectivity mesh = _state.baseMesh;
	int nverts = _state.basePos.rows();
	int nedges = _state.baseMesh.nEdges();
	int nfaces = _state.baseMesh.nFaces();

	std::set<int> pureTensionEdges, pureTensionVertices;

	getPureTensionVertsEdges(_state.tensionFaces, _state.baseMesh.faces(), &pureTensionEdges, &pureTensionVertices);

	std::set<int> clampedAmps = pureTensionVertices;

	_clampedDOFs.clear();
	if (clampedDOFs)
		_clampedDOFs = *clampedDOFs;

	std::cout << "number of clamped amp before considering pure tension case: " << _clampedDOFs.size() << std::endl;

	for (auto& vid : clampedAmps) // clampedDOFs: 0-nverts: clampedAmp, > nverts: clamped phi or dphi
	{
		if (_clampedDOFs.find(vid) == _clampedDOFs.end())
		{
			_clampedDOFs[vid] = 0;
		}
	}

	std::cout << "number of clamped amp after considering pure tension case: " << _clampedDOFs.size() << std::endl;

	std::vector<Eigen::Triplet<double> > proj;
	for (auto& eid : pureTensionEdges)
	{
		_clampedDOFs[eid + nverts] = 0;
	}

	int constrainedDOFs = _clampedDOFs.size();

	int freeDOFs = nedges + nverts - constrainedDOFs;
	int row = 0;

	for (int i = 0; i < nedges + nverts; i++)
	{
		if (_clampedDOFs.find(i) != _clampedDOFs.end())
			continue;
		if (i < nverts)
			_freeAmp++;
		else
		{
			_freeDphi++;
		}
		proj.push_back(Eigen::Triplet<double>(row, i, 1.0));
		row++;
	}
	_projM.resize(freeDOFs, nedges + nverts);
	_projM.setFromTriplets(proj.begin(), proj.end());


	std::cout << "Free Amps: " << _freeAmp << std::endl;
	std::cout << "Free dphis: " << _freeDphi << std::endl;
}


Eigen::SparseMatrix<double> TFWModel::buildIntegrabilityConstraints()
{
	int nfaces = _state.baseMesh.nFaces();
	int nedges = _state.baseMesh.nEdges();
	int nverts = _state.basePos.rows();

	Eigen::SparseMatrix<double> C;

	std::cout << "building start, faces : " << nfaces << ", edges: " << nedges << std::endl;
	std::set<int> pureTensionFaces = _state.tensionFaces;
	std::cout << "num of pure tension faces: " << pureTensionFaces.size() << std::endl;
	std::vector<Eigen::Triplet<double> > constraintCoeff;
	Eigen::SparseMatrix<double> CE, PE;
	int row = 0;
	for (int i = 0; i < nfaces; i++)
	{
		if (pureTensionFaces.find(i) != pureTensionFaces.end())
			continue;

		for (int j = 0; j < 3; j++)
		{
			int eid = _state.baseMesh.faceEdge(i, j);
			int vert1 = _state.baseMesh.edgeVertex(eid, 0);
			int vert2 = _state.baseMesh.edgeVertex(eid, 1);

			if (vert1 == _state.baseMesh.faceVertex(i, (j + 1) % 3))
				constraintCoeff.push_back(Eigen::Triplet<double>(row, eid, 1.0));
			else
				constraintCoeff.push_back(Eigen::Triplet<double>(row, eid, -1.0));
		}

		row++;
	}
	CE.resize(row, nedges);
	CE.setFromTriplets(constraintCoeff.begin(), constraintCoeff.end());
	std::cout << nfaces - row << " faces are removed from integrabilty due to pure tension check." << std::endl;
	std::cout << "Residual of integrabilty : " << (CE * _state.dphi).norm() << std::endl;

	// convert w.r.t. variable x
	PE.resize(nedges, nverts + nedges);
	PE.setZero();
	for (int i = 0; i < nedges; i++)
	{
		PE.coeffRef(i, nverts + i) = 1.0;
	}

	C = CE * PE * _projM.transpose();
	return C;

}

void TFWModel::convertCurState2Variables(const TFWState curState, Eigen::VectorXd& x)
{
	int namp = curState.amplitude.size();
	int ndphi = curState.dphi.size();
	Eigen::VectorXd y(namp + ndphi);

	//y.segment(0,namp) = curState.amplitude * _setup.rescale;
	y.segment(0, namp) = curState.amplitude;
	y.segment(namp, ndphi) = curState.dphi;

	x = _projM * y;
}

void TFWModel::convertVariables2CurState(Eigen::VectorXd x, TFWState& curState)
{
	Eigen::VectorXd fullx = _projM.transpose() * x;

	int nverts = curState.basePos.rows();
	int nedges = _state.baseMesh.nEdges();

	for (auto& it : _clampedDOFs)
	{
		fullx(it.first) = it.second;
	}

	for (int i = 0; i < nverts; i++)
	{
		curState.amplitude(i) = fullx(i);
	}
	for (int i = 0; i < nedges; i++)
	{
		curState.dphi(i) = fullx(i + nverts);
	}


}


double TFWModel::value(const Eigen::VectorXd& x)
{
	convertVariables2CurState(x, _state);
	std::shared_ptr<TFWShell> reducedShell;
	double energy = 0;
	reducedShell = std::make_shared<TFWShell>(_setup, _state, _isUsePosHess, _isParallel);
	energy = reducedShell->elasticReducedEnergy(NULL, NULL);
	return energy;

}

double TFWModel::stretchingValue(const Eigen::VectorXd& x)
{
	convertVariables2CurState(x, _state);

	// std::cout<<"initial guess: "<<_state.amplitude(0)<<std::endl;
	std::shared_ptr<TFWShell> reducedShell;
	double energy = 0;

	reducedShell = std::make_shared<TFWShell>(_setup, _state, _isUsePosHess, _isParallel);
	reducedShell->_isPosHess = _isUsePosHess;
	energy = reducedShell->stretchingEnergy(NULL, NULL);
	return energy;
}

double TFWModel::bendingValue(const Eigen::VectorXd& x)
{
	convertVariables2CurState(x, _state);

	// std::cout<<"initial guess: "<<_state.amplitude(0)<<std::endl;
	std::shared_ptr<TFWShell> reducedShell;
	double energy = 0;

	reducedShell = std::make_shared<TFWShell>(_setup, _state, _isUsePosHess, _isParallel);
	reducedShell->_isPosHess = _isUsePosHess;
	energy = reducedShell->bendingEnergy(NULL, NULL);
	return energy;
}

void TFWModel::gradient(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
{

	std::shared_ptr<TFWShell> reducedShell;
	double energy = 0;
	convertVariables2CurState(x, _state);

	reducedShell = std::make_shared<TFWShell>(_setup, _state, _isUsePosHess, _isParallel);
	reducedShell->_isPosHess = _isUsePosHess;
	energy = reducedShell->elasticReducedEnergy(&grad, NULL);

	int nverts = _state.basePos.rows();
	grad = _projM * grad;
}

void TFWModel::hessian(const Eigen::VectorXd& x, Eigen::SparseMatrix<double>& hessian)
{
	std::shared_ptr<TFWShell> reducedShell;
	double energy = 0;
	convertVariables2CurState(x, _state);
	Eigen::SparseMatrix<double> fullH;

	reducedShell = std::make_shared<TFWShell>(_setup, _state, _isUsePosHess, _isParallel);
	reducedShell->_isPosHess = _isUsePosHess;
	energy = reducedShell->elasticReducedEnergy(NULL, &fullH);


	hessian = _projM * fullH * _projM.transpose();
}

void TFWModel::save(int curIterations, TimeCost curTimeCost, double stepSize, double oldEnergy, double curEnergy, double fDelta, double ampDelta, double dphiDelta, double reg, double* projGradNorm, bool PSDHess, bool isSaveAllIntermediate)
{
	if (curIterations >= 0 && isSaveAllIntermediate)
	{
		std::string filePathPrefix = _filePrefix + "_convergence/";
		std::string energyFileName = filePathPrefix + std::string("reduced_energy.txt");
		std::ofstream efs;
		efs.open(energyFileName, std::ofstream::out | std::ofstream::app);
		bool pathExist = false;
		if (efs)
		{
			efs << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << curEnergy << std::endl;
			pathExist = true;
		}
		std::string stationResFileName = filePathPrefix + std::string("stationarity_residual.txt");
		std::ofstream srfs;
		srfs.open(stationResFileName, std::ofstream::out | std::ofstream::app);
		if (srfs && projGradNorm)
		{
			srfs << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << *projGradNorm << std::endl;
			pathExist = true;
		}
		// phi 
		std::string phiFileName = filePathPrefix + "phi/phi_" + std::to_string(curIterations) + ".txt";
		std::ofstream pfs(phiFileName);
		if (pfs)
		{
			for (int i = 0; i < _state.phi.size(); i++)
			{
				pfs << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << _state.phi(i) << std::endl;
			}
		}

		// dphi 
		std::string dphiFileName = filePathPrefix + "dphi/dphi_" + std::to_string(curIterations) + ".txt";
		std::ofstream dpfs(dphiFileName);
		if(dpfs)
		{
			for (int i = 0; i < _state.dphi.size(); i++)
			{
				dpfs << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << _state.dphi(i) << std::endl;
			}
		}
		

		// amplitude
		std::string amplitudeFileName = filePathPrefix + "amp/amp_" + std::to_string(curIterations) + ".txt";
		std::ofstream afs(amplitudeFileName);
		if(afs)
		{
			for (int i = 0; i < _state.amplitude.size(); i++)
			{
				afs << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << _state.amplitude(i) << std::endl;
			}
		}
		
		std::string ampDphiFileName = filePathPrefix + "dphi/ampWeightedDphi_" + std::to_string(curIterations) + ".csv";
		std::ofstream dapfs(ampDphiFileName);
		// amplitude weighted dphi
		if(dapfs)
		{
			for(int i = 0 ; i < _state.baseMesh.nFaces(); i++)
			{
				Eigen::Matrix2d abar = _setup.abars[i];
				Eigen::Vector3i edgeIndices;
				for (int j = 0; j < 3; j++)
				{
					edgeIndices(j) = _state.baseMesh.faceEdge(i, j);
				}

				double phiu = _state.dphi(edgeIndices(2));
				double phiv = _state.dphi(edgeIndices(1));

				int flagU = 1;
				int flagV = 1;

				if (_state.baseMesh.faceVertex(i, 0) > _state.baseMesh.faceVertex(i, 1))
				{
					flagU = -1;
				}
				if (_state.baseMesh.faceVertex(i, 0) > _state.baseMesh.faceVertex(i, 2))
				{
					flagV = -1;
				}

				phiu *= flagU;
				phiv *= flagV;
				Eigen::Vector2d w;
				w << phiu, phiv;
				w = abar.inverse() * w;

				double perfaceAmp = 0;
				for(int j = 0; j < 3; j++)
				{
					int vid = _state.baseMesh.faceVertex(i, j);
					perfaceAmp += _state.amplitude(vid) / 3.0;
				}
				w *= perfaceAmp;
				dapfs << w(0) << ",\t" << w(1) << ",\t" << 0 << ",\t" << 0 << ",\t" << 0 << ",\t" << 0 << std::endl; 
			}
		}

		// wrinkled mesh
		/*
		if(pathExist)
		{
			std::string wrinkleFileName = filePathPrefix + "wrinkledMesh/wrinkledMesh_" + std::to_string(curIterations);
			Eigen::MatrixXd NV;
			Eigen::MatrixXi NF;
			Eigen::MatrixXd upsampledTFTV, soupPhiV, soupProblemV;
			Eigen::MatrixXi upsampledTFTF, soupPhiF, soupProblemF;
			Eigen::VectorXd upsampledAmp, soupPhi;

			getUpsampledWrinkledMesh(NV, NF, upsampledTFTV, upsampledTFTF, soupPhiV, soupPhiF, soupProblemV, soupProblemF, upsampledAmp, soupPhi);
			igl::writeOBJ(wrinkleFileName + ".obj", NV, NF);

			std::string upsampledTFTpath = filePathPrefix + "upsampled/upsampledTFT_" + std::to_string(curIterations) + ".obj";
			igl::writeOBJ(upsampledTFTpath, upsampledTFTV, upsampledTFTF);

			std::string soupPath = filePathPrefix + "upsampled/phiSoup_" + std::to_string(curIterations) + ".obj";
			igl::writeOBJ(soupPath, soupPhiV, soupPhiF);

			std::string upsampledAmpPath = filePathPrefix + "upsampled/upsampledAmp_" + std::to_string(curIterations) + ".csv";
			std::ofstream upsampledAmpFile(upsampledAmpPath);
			for (int i = 0; i < upsampledAmp.size(); i++)
				upsampledAmpFile << upsampledAmp[i] << ",\t" << 3.14159 << std::endl; // fix for buggy Houdini import

			std::string upsampledPhiPath = filePathPrefix + "upsampled/phiSoup_" + std::to_string(curIterations) + ".csv";
			std::ofstream upsampledPhiFile(upsampledPhiPath);
			for (int i = 0; i < soupPhi.size(); i++)
				upsampledPhiFile << soupPhi[i] << ",\t" << 3.14159 << std::endl;

			std::string problempath = filePathPrefix + "upsampled/problemV_" + std::to_string(curIterations) + ".obj";
			igl::writeOBJ(problempath, soupProblemV, soupProblemF);

		}
		*/
	}

	if (curIterations == 0)
		return; // initial state
	std::string statusFileName = _filePrefix + std::string("_reduced_status.txt");
	std::ofstream sfs;
	sfs.open(statusFileName, std::ofstream::out | std::ofstream::app);

	if (sfs)
	{
		sfs << std::endl << "iter: " << curIterations << ", total time: " << curTimeCost.totalTime() << std::endl;
		sfs << "line search rate: " << stepSize << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << ", actual hessian: "; 
		if (PSDHess)
			sfs << "false" << ", reg: " << reg << std::endl;
		else
			sfs << "true" << ", reg: " << reg << std::endl;

		sfs << "f_old: " << oldEnergy << ", f_new: " << curEnergy << ", f_delta: " << fDelta << ", amp_delta: " << ampDelta << ", dphi_delta: " << dphiDelta;
		if (projGradNorm)
			sfs <<", ||stationarity_residual|| = " <<*(projGradNorm) << std::endl;
		else
			sfs << std::endl;
	}


	std::string timeFileName = _filePrefix + std::string("_reduced_timing.txt");
	std::ofstream tfs;
	tfs.open(timeFileName, std::ofstream::out | std::ofstream::app);

	if (tfs)
	{
		tfs << "iter: " << curIterations << ", total time: " << curTimeCost.totalTime() << std::endl;
		tfs << "total grad_cost: " << curTimeCost.gradTime << ", total hess_cost: " << curTimeCost.hessTime << ", total solver_cost: " << curTimeCost.solverTime << ", total LLT_cost: " << curTimeCost.unconstrainedLLTTime << ", total linesearch_cost: " << curTimeCost.lineSearchTime << ", total update_cost: " << curTimeCost.updateTime << ", total checkConverg_cost: " << curTimeCost.convergenceCheckTime << std::endl;
	}

	int num = 10;
	if (curIterations % num == 0)
	{
		std::stringstream Filename;
		Filename << _filePrefix << "_cur";
		std::cout<<"save file in "<<Filename.str()<<std::endl;
		std::string prefix = Filename.str();

		// phi 
		std::string phiFileName = prefix + std::string("_phi.txt");
		std::ofstream pfs(phiFileName);
		for(int i = 0; i < _state.phi.size(); i++)
		{
			pfs << std::setprecision(std::numeric_limits<long double>::digits10 + 1)<<_state.phi(i)<<std::endl;
		}

		// dphi 
		std::string dphiFileName = prefix + std::string("_dphi.txt");
		std::ofstream dpfs(dphiFileName);
		for(int i = 0; i < _state.dphi.size(); i++)
		{
			dpfs << std::setprecision(std::numeric_limits<long double>::digits10 + 1)<<_state.dphi(i)<<std::endl;
		}

		// amplitude
		std::string amplitudeFileName = prefix + std::string("_amplitude.txt");
		std::ofstream afs(amplitudeFileName);
		afs.precision(15);
		for(int i = 0; i < _state.amplitude.size(); i++)
		{
			afs <<std::setprecision(std::numeric_limits<long double>::digits10 + 1)<<_state.amplitude(i)<<std::endl;
		}
	}
}

void TFWModel::getUpsampledWrinkledMesh(Eigen::MatrixXd& NV, Eigen::MatrixXi& NF, Eigen::MatrixXd &upsampledTFTV, Eigen::MatrixXi &upsampledTFTF, Eigen::MatrixXd &soupPhiV, Eigen::MatrixXi &soupPhiF, Eigen::MatrixXd &soupProblemV, Eigen::MatrixXi &soupProblemF, Eigen::VectorXd &upsampledAmp, Eigen::VectorXd &soupPhi)
{
	auto curState = _state;
	Eigen::MatrixXd curV = curState.basePos;
	Eigen::MatrixXi curF = curState.baseMesh.faces();
	Eigen::VectorXd curAmp = curState.amplitude;
	Eigen::VectorXd curPhi = curState.phi;

	std::set<int> clampedVerts;

	for (auto& it : _setup.clampedDOFs)
	{
		int vid = it.first / 3;
		if (clampedVerts.count(vid) == 0)
			clampedVerts.insert(vid);
	}
	std::set<int> problemFaces;
	auto sff = std::make_shared<MidedgeAverageFormulation>();

	roundPhiFromDphiCutbyTension(curState.basePos, curState.baseMesh.faces(), curV, curF, _setup.abars, curState.amplitude, curState.dphi, GurobiRound, curState.phi, curPhi, curAmp, problemFaces);
	wrinkledMeshUpsamplingUncut(curState.basePos, curState.baseMesh.faces(), _setup.restV, _setup.restF, curV, curF, problemFaces, clampedVerts,
		&NV, &NF, &upsampledTFTV, &upsampledTFTF, &soupPhiV, &soupPhiF, &soupProblemV, &soupProblemF, &upsampledAmp, &soupPhi,
		curAmp, curPhi, *sff, _setup.YoungsModulus, _setup.PoissonsRatio, 2, SubdivisionType::Loop, false, true);
}

Eigen::VectorXd TFWModel::getProjectedGradient(const Eigen::VectorXd &x)
{
	Eigen::VectorXd g;
	Eigen::VectorXd Beq, BIneq;
	Eigen::SparseMatrix<double> Aeq, AIneq, I;
	Eigen::VectorXd lx, ux;

	gradient(x, g);
	return g;

//	I.resize(g.rows(), g.rows());
//	I.setIdentity();
//
//	Aeq = buildIntegrabilityConstraints();
//	Beq = -Aeq * x;
//	Beq.setZero();
//
//	AIneq.resize(0, 0);
//	BIneq.resize(0);
//
//	lx = -x;
//	// for(int i = 0; i < _freeAmp; i++)
//	//     lx(i) = x(i) > 0 ? -x(i) : 0;
//	lx.segment(_freeAmp, x.size() - _freeAmp).setConstant(-std::numeric_limits<double>::infinity());
//	ux.resize(0);
//
//	EigenNASOQSparse solver;
//	solver.setAccThresh(1e-10);
//
//	Eigen::VectorXd B = g;
//
//	Eigen::VectorXd projGrad;
//	double perturb = 1e-9;
//	solver.solve(I, B, Aeq, Beq, AIneq, BIneq, lx, ux, projGrad, perturb);
//
//	return projGrad;
}

void TFWModel::testValueAndGradient(const Eigen::VectorXd &x)
{
	std::cout << "Test value and gradient. " << std::endl;
	double f = value(x);
	Eigen::VectorXd grad;
	gradient(x, grad);
	Eigen::VectorXd dir = Eigen::VectorXd::Random(x.size());
	dir.normalize();
//    Eigen::VectorXd finiteDiff;
//    finiteGradient(x, finiteDiff);
//    std::cout<<grad - finiteDiff<<std::endl;
//    std::cout<<"Error norm: "<<(grad - finiteDiff).lpNorm<Eigen::Infinity>()<<std::endl;
//    std::cout<<std::endl<<"grad \t difference " <<std::endl;
//    for (int i = 0; i < grad.size(); i++)
//        std::cout << grad(i) << " " << grad(i) - finiteDiff(i) << std::endl;
	for(int i = 3; i < 10; i++)
	{
		double eps = std::pow(10, -i);
		Eigen::VectorXd x1 = x + eps * dir;
		Eigen::VectorXd x2 = x - eps * dir;
		double f1 = value(x1);
		double f2 = value(x2);
		std::cout<<std::endl<<"eps: "<<eps<<std::endl;
		std::cout<<std::setprecision(std::numeric_limits<long double>::digits10 + 1)<<"energy: "<<f<<", energy after perturbation (right, left): "<<f1<<", "<<f2<<std::endl;
		std::cout<<std::setprecision(6);
		std::cout<<"right finite difference: "<<(f1 - f) / eps<<std::endl;
		std::cout<<"left finite difference: "<<(f - f2) / eps <<std::endl;
		std::cout<<"central difference: "<<(f1 - f2) / 2 / eps<<std::endl;
		std::cout<<"direction derivative: "<<grad.dot(dir)<<std::endl;
		std::cout<<"right error: "<<std::abs((f1 - f) / eps - grad.dot(dir))<<", left error: "<<std::abs((f - f2) / eps - grad.dot(dir))<<", central error: "<<std::abs((f1 - f2) / 2 / eps - grad.dot(dir))<<std::endl;
	}
}

void TFWModel::testGradientAndHessian(const Eigen::VectorXd& x)
{
	bool isUsePosHess = _isUsePosHess;
	_isUsePosHess = false;
	std::cout << "Test gradient and hessian. " << std::endl;
	Eigen::SparseMatrix<double> hess;
	Eigen::VectorXd deriv;
	gradient(x, deriv);
	hessian(x, hess);

	Eigen::VectorXd dir = Eigen::VectorXd::Random(x.size());
	dir.normalize();
	for (int i = 3; i < 10; i++)
	{
		double eps = std::pow(10, -i);
		Eigen::VectorXd x1 = x + eps * dir;
		Eigen::VectorXd deriv1;
		gradient(x1, deriv1);

		std::cout << std::endl << "eps: " << eps << std::endl;
		std::cout << std::setprecision(6);
		std::cout << "finite difference: " << (deriv1 - deriv).norm() / eps << std::endl;
		std::cout << "direction derivative: " << (hess * dir).norm() << std::endl;
		std::cout << "error: " << ((deriv1 - deriv) / eps - hess * dir).norm() << std::endl;
	}
	_isUsePosHess = isUsePosHess;
}


void TFWModel::checkPD4ConstraintedHess(const Eigen::VectorXd& x)
{
	bool isUsePosHess = _isUsePosHess;
	_isUsePosHess = false;
	std::cout << "Test constriant H PD. " << std::endl;
	Eigen::SparseMatrix<double> hess;
	Eigen::VectorXd deriv;
	gradient(x, deriv);
	hessian(x, hess);

	std::cout << "check PD for entire hessian" << std::endl;
	Eigen::SparseMatrix<double> idmat(hess.rows(), hess.cols());
	idmat.setIdentity();
	Eigen::SparseMatrix<double> H = hess;

	Eigen::CholmodSimplicialLLT<Eigen::SparseMatrix<double> > solver0(H);
	if (solver0.info() == Eigen::Success)
	{
		std::cout << "entire hessian is PD." << std::endl;
		return;
	}
	double reg = 1e-10;

	while (solver0.info() != Eigen::Success)
	{
		reg = std::max(2 * reg, 1e-16);
		std::cout << "Matrix is not positive definite, current reg = " << reg << std::endl;
		H = hess + reg * idmat;
		solver0.compute(H);
	}


	Eigen::SparseMatrix<double> C = buildIntegrabilityConstraints();    // C * x = 0
	Eigen::MatrixXd N = nullspaceExtraction(C);

	std::cout << "Null space check: " << (C * N).norm() << std::endl;


	Eigen::SparseMatrix<double> constraintedH = (N.transpose() * hess * N).sparseView();


	Eigen::SparseMatrix<double> constraintedH_reg = constraintedH;
	Eigen::SparseMatrix<double> I(constraintedH.rows(), constraintedH.cols());
	I.setIdentity();

	Eigen::CholmodSimplicialLLT<Eigen::SparseMatrix<double> > solver(constraintedH_reg);
	if (solver.info() == Eigen::Success)
	{
		std::cout << "Constraint hessian is PD." << std::endl;
		return;
	}
	reg = 1e-10;

	while (solver.info() != Eigen::Success)
	{
		reg = std::max(2 * reg, 1e-16);
		std::cout << "Matrix is not positive definite, current reg = " << reg << std::endl;
		constraintedH_reg = constraintedH + reg * I;
		solver.compute(constraintedH_reg);
	}

	_isUsePosHess = isUsePosHess;
}