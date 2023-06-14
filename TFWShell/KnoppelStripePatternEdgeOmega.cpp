#include "KnoppelStripePatternEdgeOmega.h"
#include <igl/cotmatrix.h>
#include <SymGEigsShiftSolver.h>
#include <MatOp/SparseCholesky.h>
#include <Eigen/CholmodSupport>
#include <MatOp/SparseSymShiftSolve.h>
#include <iostream>

void computeEdgeMatrix(const MeshConnectivity &mesh, const Eigen::VectorXd& edgeW, const Eigen::VectorXd& edgeWeight,
									 const int nverts, Eigen::SparseMatrix<double> &A)
{
	std::vector<Eigen::Triplet<double>> AT;
	int nfaces = mesh.nFaces();
	int nedges = mesh.nEdges();

	for(int i = 0; i < nedges; i++)
	{
		int vid0 = mesh.edgeVertex(i, 0);
		int vid1 = mesh.edgeVertex(i, 1);

		AT.push_back({2 * vid0, 2 * vid0, edgeWeight(i)});
		AT.push_back({2 * vid0 + 1, 2 * vid0 + 1, edgeWeight(i)});

		AT.push_back({2 * vid1, 2 * vid1, edgeWeight(i)});
		AT.push_back({2 * vid1 + 1, 2 * vid1 + 1, edgeWeight(i)});

		std::complex<double> expw0 = std::complex<double>(std::cos(edgeW(i)), std::sin(edgeW(i)));

		AT.push_back({2 * vid0, 2 * vid1, -edgeWeight(i) * (expw0.real())});
		AT.push_back({2 * vid0 + 1, 2 * vid1, -edgeWeight(i) * (-expw0.imag())});
		AT.push_back({2 * vid0, 2 * vid1 + 1, -edgeWeight(i) * (expw0.imag())});
		AT.push_back({2 * vid0 + 1, 2 * vid1 + 1, -edgeWeight(i) * (expw0.real())});

		AT.push_back({ 2 * vid1, 2 * vid0, -edgeWeight(i) * (expw0.real()) });
		AT.push_back({ 2 * vid1, 2 * vid0 + 1, -edgeWeight(i) * (-expw0.imag()) });
		AT.push_back({ 2 * vid1 + 1, 2 * vid0, -edgeWeight(i) * (expw0.imag()) });
		AT.push_back({ 2 * vid1 + 1, 2 * vid0 + 1, -edgeWeight(i) * (expw0.real()) });

	}
	A.resize(2 * nverts, 2 * nverts);
	A.setFromTriplets(AT.begin(), AT.end());
}

void roundZvalsFromEdgeOmega(const MeshConnectivity &mesh, const Eigen::VectorXd& edgeW,
	const Eigen::VectorXd& edgeWeight, const Eigen::VectorXd& vertArea, const int nverts, std::vector<std::complex<double>> &zvals)
{
	std::vector<Eigen::Triplet<double>> BT;
	int nfaces = mesh.nFaces();
	int nedges = mesh.nEdges();

	for (int i = 0; i < nverts; i++)
	{
		BT.push_back({ 2 * i, 2 * i, vertArea(i) });
		BT.push_back({ 2 * i + 1, 2 * i + 1, vertArea(i) });
	}
	
	
	Eigen::SparseMatrix<double> A;
	computeEdgeMatrix(mesh, edgeW, edgeWeight, nverts, A);

	Eigen::SparseMatrix<double> B(2 * nverts, 2 * nverts);
	B.setFromTriplets(BT.begin(), BT.end());
   /* std::cout << A.toDense() << std::endl;
	std::cout << B.toDense() << std::endl;*/

	Spectra::SymShiftInvert<double> op(A, B);
	Spectra::SparseSymMatProd<double> Bop(B);
	Spectra::SymGEigsShiftSolver<Spectra::SymShiftInvert<double>, Spectra::SparseSymMatProd<double>, Spectra::GEigsMode::ShiftInvert> geigs(op, Bop, 1, 6, -1e-6);
	geigs.init();
	int nconv = geigs.compute(Spectra::SortRule::LargestMagn, 1e6);

	Eigen::VectorXd evalues;
	Eigen::MatrixXd evecs;

	evalues = geigs.eigenvalues();
	evecs = geigs.eigenvectors();
	if (nconv != 1 || geigs.info() != Spectra::CompInfo::Successful)
	{
		std::cout << "Eigensolver failed to converge!!" << std::endl;
	}

	std::cout << "Eigenvalue is " << evalues[0] << std::endl;

	zvals.clear();
	for(int i = 0; i < nverts; i++)
	{
		zvals.push_back(std::complex<double>(evecs(2 * i, 0), evecs(2 * i + 1, 0)));
	}
}


