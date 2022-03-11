
#include <igl/cotmatrix.h>
#include <igl/boundary_facets.h>
#include <igl/boundary_loop.h>
#include "WTFShell.h"

using namespace WTF;

double WTFShell::computeAmplitudesFromQuad(int faceId, int quadId, Eigen::Vector2d* da, Eigen::Vector3d* gradA, Eigen::Matrix<double, 2, 3>* gradDA, Eigen::Matrix<double, 3, 3>* hessianA, std::vector<Eigen::Matrix<double, 3, 3>>* hessianDA)
{
	Eigen::Vector3d amps;
	Eigen::Vector3i vertices;
	for (int i = 0; i < 3; i++)
	{
		vertices(i) = _state.baseMesh.faceVertex(faceId, i);
		amps(i) = _state.amplitude(vertices(i));
	}
	double u = _setup.quadPoints[quadId].u;
	double v = _setup.quadPoints[quadId].v;
	double a = u * amps(1) + v * amps(2) + (1 - u - v) * amps(0);


	if (da)
	{
		da->coeffRef(0) = amps(1) - amps(0);
		da->coeffRef(1) = amps(2) - amps(0);
	}

	if (gradA)
	{
		gradA->coeffRef(0) = 1 - u - v;
		gradA->coeffRef(1) = u;
		gradA->coeffRef(2) = v;
	}

	if (hessianA)
	{
		hessianA->setZero();
	}

	if (gradDA)
	{
		gradDA->coeffRef(0, 0) = -1.0;
		gradDA->coeffRef(0, 1) = 1.0;
		gradDA->coeffRef(0, 2) = 0;

		gradDA->coeffRef(1, 0) = -1.0;
		gradDA->coeffRef(1, 1) = 0;
		gradDA->coeffRef(1, 2) = 1.0;

	}
	if (hessianDA)
	{
		hessianDA->resize(2);
		hessianDA->at(0).setZero();
		hessianDA->at(1).setZero();
	}

	return a;
}

Eigen::Vector2d WTFShell::computeDphi(int faceId, Eigen::Matrix<double, 2, 3>* gradDphi)
{
	Eigen::Vector3i edgeIndices;
	for (int i = 0; i < 3; i++)
	{
		edgeIndices(i) = _state.baseMesh.faceEdge(faceId, i);
	}

	double phiu = _state.dphi(edgeIndices(2));
	double phiv = _state.dphi(edgeIndices(1));

	int flagU = 1;
	int flagV = 1;

	if (_state.baseMesh.faceVertex(faceId, 0) > _state.baseMesh.faceVertex(faceId, 1))
	{
		flagU = -1;
	}
	if (_state.baseMesh.faceVertex(faceId, 0) > _state.baseMesh.faceVertex(faceId, 2))
	{
		flagV = -1;
	}

	phiu *= flagU;
	phiv *= flagV;

	if (gradDphi)
	{
		gradDphi->setZero();
		gradDphi->coeffRef(0, 2) = flagU;
		gradDphi->coeffRef(1, 1) = flagV;
	}

	return Eigen::Vector2d(phiu, phiv);
}

Eigen::Matrix2d WTFShell::computeDaDphiTensor(int faceId, int quadId, std::vector<Eigen::Matrix2d>* deriv, std::vector<Eigen::Matrix<double, 6, 6>>* hessian)
{
	Eigen::Vector2d da, dphi;
	Eigen::Vector3d gradA;
	Eigen::Matrix<double, 2, 3> gradDA, gradDphi;
	Eigen::Matrix<double, 3, 3> hessA;
	std::vector<Eigen::Matrix<double, 3, 3> > hessDA;

	double a = computeAmplitudesFromQuad(faceId, quadId, &da, NULL, deriv ? &gradDA : NULL, hessian ? &hessA : NULL, hessian ? &hessDA : NULL);
	dphi = computeDphi(faceId, deriv ? &gradDphi : NULL);

	if (deriv)
	{
		for (int i = 0; i < 3; i++)
		{
			Eigen::Matrix2d grad_dadphiT;
			grad_dadphiT(0, 0) = gradDA(0, i) * dphi(0);
			grad_dadphiT(0, 1) = gradDA(0, i) * dphi(1);
			grad_dadphiT(1, 0) = gradDA(1, i) * dphi(0);
			grad_dadphiT(1, 1) = gradDA(1, i) * dphi(1);

			deriv->push_back(grad_dadphiT);
		}

		for (int i = 0; i < 3; i++)
		{
			Eigen::Matrix2d grad_dadphiT;
			grad_dadphiT(0, 0) = gradDphi(0, i) * da(0);
			grad_dadphiT(0, 1) = gradDphi(1, i) * da(0);
			grad_dadphiT(1, 0) = gradDphi(0, i) * da(1);
			grad_dadphiT(1, 1) = gradDphi(1, i) * da(1);

			deriv->push_back(grad_dadphiT);
		}
	}

	if (hessian)
	{
		hessian->resize(4);

		for (int i = 0; i < 4; i++)
			hessian->at(i).setZero();

		hessian->at(0).block(0, 0, 3, 3) << hessDA[0](0, 0) * dphi(0), hessDA[0](0, 1)* dphi(0), hessDA[0](0, 2)* dphi(0),
			hessDA[0](1, 0)* dphi(0), hessDA[0](1, 1)* dphi(0), hessDA[0](1, 2)* dphi(0),
			hessDA[0](2, 0)* dphi(0), hessDA[0](2, 1)* dphi(0), hessDA[0](2, 2)* dphi(0);
		hessian->at(0).block(0, 3, 3, 3) <<
			gradDA(0, 0) * gradDphi(0, 0), gradDA(0, 0)* gradDphi(0, 1), gradDA(0, 0)* gradDphi(0, 2),
			gradDA(0, 1)* gradDphi(0, 0), gradDA(0, 1)* gradDphi(0, 1), gradDA(0, 1)* gradDphi(0, 2),
			gradDA(0, 2)* gradDphi(0, 0), gradDA(0, 2)* gradDphi(0, 1), gradDA(0, 2)* gradDphi(0, 2);
		hessian->at(0).block(3, 0, 3, 3) = hessian->at(0).block(0, 3, 3, 3).transpose();


		hessian->at(1).block(0, 0, 3, 3) << hessDA[0](0, 0) * dphi(1), hessDA[0](0, 1)* dphi(1), hessDA[0](0, 2)* dphi(1),
			hessDA[0](1, 0)* dphi(1), hessDA[0](1, 1)* dphi(1), hessDA[0](1, 2)* dphi(1),
			hessDA[0](2, 0)* dphi(1), hessDA[0](2, 1)* dphi(1), hessDA[0](2, 2)* dphi(1);
		hessian->at(1).block(0, 3, 3, 3) <<
			gradDA(0, 0) * gradDphi(1, 0), gradDA(0, 0)* gradDphi(1, 1), gradDA(0, 0)* gradDphi(1, 2),
			gradDA(0, 1)* gradDphi(1, 0), gradDA(0, 1)* gradDphi(1, 1), gradDA(0, 1)* gradDphi(1, 2),
			gradDA(0, 2)* gradDphi(1, 0), gradDA(0, 2)* gradDphi(1, 1), gradDA(0, 2)* gradDphi(1, 2);
		hessian->at(1).block(3, 0, 3, 3) = hessian->at(1).block(0, 3, 3, 3).transpose();


		hessian->at(2).block(0, 0, 3, 3) << hessDA[1](0, 0) * dphi(0), hessDA[1](0, 1)* dphi(0), hessDA[1](0, 2)* dphi(0),
			hessDA[1](1, 0)* dphi(0), hessDA[1](1, 1)* dphi(0), hessDA[1](1, 2)* dphi(0),
			hessDA[1](2, 0)* dphi(0), hessDA[1](2, 1)* dphi(0), hessDA[1](2, 2)* dphi(0);
		hessian->at(2).block(0, 3, 3, 3) <<
			gradDA(1, 0) * gradDphi(0, 0), gradDA(1, 0)* gradDphi(0, 1), gradDA(1, 0)* gradDphi(0, 2),
			gradDA(1, 1)* gradDphi(0, 0), gradDA(1, 1)* gradDphi(0, 1), gradDA(1, 1)* gradDphi(0, 2),
			gradDA(1, 2)* gradDphi(0, 0), gradDA(1, 2)* gradDphi(0, 1), gradDA(1, 2)* gradDphi(0, 2);
		hessian->at(2).block(3, 0, 3, 3) = hessian->at(2).block(0, 3, 3, 3).transpose();


		hessian->at(3).block(0, 0, 3, 3) << hessDA[1](0, 0) * dphi(1), hessDA[1](0, 1)* dphi(1), hessDA[1](0, 2)* dphi(1),
			hessDA[1](1, 0)* dphi(1), hessDA[1](1, 1)* dphi(1), hessDA[1](1, 2)* dphi(1),
			hessDA[1](2, 0)* dphi(1), hessDA[1](2, 1)* dphi(1), hessDA[1](2, 2)* dphi(1);
		hessian->at(3).block(0, 3, 3, 3) <<
			gradDA(1, 0) * gradDphi(1, 0), gradDA(1, 0)* gradDphi(1, 1), gradDA(1, 0)* gradDphi(1, 2),
			gradDA(1, 1)* gradDphi(1, 0), gradDA(1, 1)* gradDphi(1, 1), gradDA(1, 1)* gradDphi(1, 2),
			gradDA(1, 2)* gradDphi(1, 0), gradDA(1, 2)* gradDphi(1, 1), gradDA(1, 2)* gradDphi(1, 2);
		hessian->at(3).block(3, 0, 3, 3) = hessian->at(3).block(0, 3, 3, 3).transpose();

	}

	return da * dphi.transpose();
}

Eigen::Matrix2d WTFShell::computeDphiDaTensor(int faceId, int quadId, std::vector<Eigen::Matrix2d>* deriv, std::vector<Eigen::Matrix<double, 6, 6>>* hessian)
{
	std::vector<Eigen::Matrix2d> deriv1;
	std::vector<Eigen::Matrix<double, 6, 6>> hessian1;

	Eigen::Matrix2d dadphTensor = computeDaDphiTensor(faceId, quadId, deriv ? &deriv1 : NULL, hessian ? &hessian1 : NULL);

	if (deriv)
	{
		for (int i = 0; i < deriv1.size(); i++)
			deriv->push_back(deriv1[i].transpose());
	}

	if (hessian)
	{
		hessian->resize(4);
		hessian->at(0) = hessian1[0];
		hessian->at(1) = hessian1[2];
		hessian->at(2) = hessian1[1];
		hessian->at(3) = hessian1[3];
	}

	return dadphTensor.transpose();
}

Eigen::Matrix2d WTFShell::computeDaDaTensor(int faceId, int quadId, std::vector<Eigen::Matrix2d>* deriv, std::vector<Eigen::Matrix<double, 3, 3>>* hessian)
{
	Eigen::Vector2d da;
	Eigen::Vector3d gradA;
	Eigen::Matrix<double, 2, 3> gradDA;
	Eigen::Matrix<double, 3, 3> hessA;
	std::vector<Eigen::Matrix<double, 3, 3> > hessDA;

	double a = computeAmplitudesFromQuad(faceId, quadId, &da, NULL, deriv ? &gradDA : NULL, hessian ? &hessA : NULL, hessian ? &hessDA : NULL);

	if (deriv)
	{
		for (int i = 0; i < 3; i++)
		{
			Eigen::Matrix2d grad_dadaT;
			grad_dadaT(0, 0) = 2 * da(0) * gradDA(0, i);
			grad_dadaT(0, 1) = da(0) * gradDA(1, i) + da(1) * gradDA(0, i);
			grad_dadaT(1, 0) = grad_dadaT(0, 1);
			grad_dadaT(1, 1) = 2 * da(1) * gradDA(1, i);

			deriv->push_back(grad_dadaT);
		}
	}

	if (hessian)
	{
		hessian->resize(4);
		hessian->at(0) <<
			2 * gradDA(0, 0) * gradDA(0, 0) + 2.0 * da(0) * hessDA[0](0, 0), 2 * gradDA(0, 1) * gradDA(0, 0) + 2.0 * da(0) * hessDA[0](0, 1), 2 * gradDA(0, 2) * gradDA(0, 0) + 2.0 * da(0) * hessDA[0](0, 2),
			2 * gradDA(0, 0) * gradDA(0, 1) + 2.0 * da(0) * hessDA[0](1, 0), 2 * gradDA(0, 1) * gradDA(0, 1) + 2.0 * da(0) * hessDA[0](1, 1), 2 * gradDA(0, 2) * gradDA(0, 1) + 2.0 * da(0) * hessDA[0](1, 2),
			2 * gradDA(0, 0) * gradDA(0, 2) + 2.0 * da(0) * hessDA[0](2, 0), 2 * gradDA(0, 1) * gradDA(0, 2) + 2.0 * da(0) * hessDA[0](2, 1), 2 * gradDA(0, 2) * gradDA(0, 2) + 2.0 * da(0) * hessDA[0](2, 2);

		hessian->at(1) <<
			gradDA(0, 0) * gradDA(1, 0) + gradDA(1, 0) * gradDA(0, 0) + da(0) * hessDA[1](0, 0) + da(1) * hessDA[0](0, 0), gradDA(0, 1)* gradDA(1, 0) + gradDA(1, 1) * gradDA(0, 0) + da(0) * hessDA[1](0, 1) + da(1) * hessDA[0](0, 1), gradDA(0, 2)* gradDA(1, 0) + gradDA(1, 2) * gradDA(0, 0) + da(0) * hessDA[1](0, 2) + da(1) * hessDA[0](0, 2),
			gradDA(0, 0)* gradDA(1, 1) + gradDA(1, 0) * gradDA(0, 1) + da(0) * hessDA[1](1, 0) + da(1) * hessDA[0](1, 0), gradDA(0, 1)* gradDA(1, 1) + gradDA(1, 1) * gradDA(0, 1) + da(0) * hessDA[1](1, 1) + da(1) * hessDA[0](1, 1), gradDA(0, 2)* gradDA(1, 1) + gradDA(1, 2) * gradDA(0, 1) + da(0) * hessDA[1](1, 2) + da(1) * hessDA[0](1, 2),
			gradDA(0, 0)* gradDA(1, 2) + gradDA(1, 0) * gradDA(0, 2) + da(0) * hessDA[1](2, 0) + da(1) * hessDA[0](2, 0), gradDA(0, 1)* gradDA(1, 2) + gradDA(1, 1) * gradDA(0, 2) + da(0) * hessDA[1](2, 1) + da(1) * hessDA[0](2, 1), gradDA(0, 2)* gradDA(1, 2) + gradDA(1, 2) * gradDA(0, 2) + da(0) * hessDA[1](2, 2) + da(1) * hessDA[0](2, 2);

		hessian->at(2) <<
			gradDA(0, 0) * gradDA(1, 0) + gradDA(1, 0) * gradDA(0, 0) + da(0) * hessDA[1](0, 0) + da(1) * hessDA[0](0, 0), gradDA(0, 1)* gradDA(1, 0) + gradDA(1, 1) * gradDA(0, 0) + da(0) * hessDA[1](0, 1) + da(1) * hessDA[0](0, 1), gradDA(0, 2)* gradDA(1, 0) + gradDA(1, 2) * gradDA(0, 0) + da(0) * hessDA[1](0, 2) + da(1) * hessDA[0](0, 2),
			gradDA(0, 0)* gradDA(1, 1) + gradDA(1, 0) * gradDA(0, 1) + da(0) * hessDA[1](1, 0) + da(1) * hessDA[0](1, 0), gradDA(0, 1)* gradDA(1, 1) + gradDA(1, 1) * gradDA(0, 1) + da(0) * hessDA[1](1, 1) + da(1) * hessDA[0](1, 1), gradDA(0, 2)* gradDA(1, 1) + gradDA(1, 2) * gradDA(0, 1) + da(0) * hessDA[1](1, 2) + da(1) * hessDA[0](1, 2),
			gradDA(0, 0)* gradDA(1, 2) + gradDA(1, 0) * gradDA(0, 2) + da(0) * hessDA[1](2, 0) + da(1) * hessDA[0](2, 0), gradDA(0, 1)* gradDA(1, 2) + gradDA(1, 1) * gradDA(0, 2) + da(0) * hessDA[1](2, 1) + da(1) * hessDA[0](2, 1), gradDA(0, 2)* gradDA(1, 2) + gradDA(1, 2) * gradDA(0, 2) + da(0) * hessDA[1](2, 2) + da(1) * hessDA[0](2, 2);

		hessian->at(3) <<
			2 * gradDA(1, 0) * gradDA(1, 0) + 2.0 * da(1) * hessDA[1](0, 0), 2 * gradDA(1, 1) * gradDA(1, 0) + 2.0 * da(1) * hessDA[1](0, 1), 2 * gradDA(1, 2) * gradDA(1, 0) + 2.0 * da(1) * hessDA[1](0, 2),
			2 * gradDA(1, 0) * gradDA(1, 1) + 2.0 * da(1) * hessDA[1](1, 0), 2 * gradDA(1, 1) * gradDA(1, 1) + 2.0 * da(1) * hessDA[1](1, 1), 2 * gradDA(1, 2) * gradDA(1, 1) + 2.0 * da(1) * hessDA[1](1, 2),
			2 * gradDA(1, 0) * gradDA(1, 2) + 2.0 * da(1) * hessDA[1](2, 0), 2 * gradDA(1, 1) * gradDA(1, 2) + 2.0 * da(1) * hessDA[1](2, 1), 2 * gradDA(1, 2) * gradDA(1, 2) + 2.0 * da(1) * hessDA[1](2, 2);
	}

	return da * da.transpose();
}

Eigen::Matrix2d WTFShell::computeDphiDphiTensor(int faceId, int quadId, std::vector<Eigen::Matrix2d>* deriv, std::vector<Eigen::Matrix<double, 3, 3>>* hessian)
{
	Eigen::Vector2d dphi;
	Eigen::Matrix<double, 2, 3> gradDphi;
	dphi = computeDphi(faceId, deriv ? &gradDphi : NULL);

	if (deriv)
	{
		for (int i = 0; i < 3; i++)
		{
			Eigen::Matrix2d grad_dphidphiT;
			grad_dphidphiT(0, 0) = 2 * dphi(0) * gradDphi(0, i);
			grad_dphidphiT(0, 1) = dphi(0) * gradDphi(1, i) + dphi(1) * gradDphi(0, i);
			grad_dphidphiT(1, 0) = grad_dphidphiT(0, 1);
			grad_dphidphiT(1, 1) = 2 * dphi(1) * gradDphi(1, i);

			deriv->push_back(grad_dphidphiT);
		}
	}

	if (hessian)
	{
		hessian->resize(4);
		hessian->at(0) <<
			2 * gradDphi(0, 0) * gradDphi(0, 0), 2 * gradDphi(0, 1) * gradDphi(0, 0), 2 * gradDphi(0, 2) * gradDphi(0, 0),
			2 * gradDphi(0, 0) * gradDphi(0, 1), 2 * gradDphi(0, 1) * gradDphi(0, 1), 2 * gradDphi(0, 2) * gradDphi(0, 1),
			2 * gradDphi(0, 0) * gradDphi(0, 2), 2 * gradDphi(0, 1) * gradDphi(0, 2), 2 * gradDphi(0, 2) * gradDphi(0, 2);

		hessian->at(1) <<
			gradDphi(0, 0) * gradDphi(1, 0) + gradDphi(1, 0) * gradDphi(0, 0), gradDphi(0, 1)* gradDphi(1, 0) + gradDphi(1, 1) * gradDphi(0, 0), gradDphi(0, 2)* gradDphi(1, 0) + gradDphi(1, 2) * gradDphi(0, 0),
			gradDphi(0, 0)* gradDphi(1, 1) + gradDphi(1, 0) * gradDphi(0, 1), gradDphi(0, 1)* gradDphi(1, 1) + gradDphi(1, 1) * gradDphi(0, 1), gradDphi(0, 2)* gradDphi(1, 1) + gradDphi(1, 2) * gradDphi(0, 1),
			gradDphi(0, 0)* gradDphi(1, 2) + gradDphi(1, 0) * gradDphi(0, 2), gradDphi(0, 1)* gradDphi(1, 2) + gradDphi(1, 1) * gradDphi(0, 2), gradDphi(0, 2)* gradDphi(1, 2) + gradDphi(1, 2) * gradDphi(0, 2);

		hessian->at(2) <<
			gradDphi(0, 0) * gradDphi(1, 0) + gradDphi(1, 0) * gradDphi(0, 0), gradDphi(0, 1)* gradDphi(1, 0) + gradDphi(1, 1) * gradDphi(0, 0), gradDphi(0, 2)* gradDphi(1, 0) + gradDphi(1, 2) * gradDphi(0, 0),
			gradDphi(0, 0)* gradDphi(1, 1) + gradDphi(1, 0) * gradDphi(0, 1), gradDphi(0, 1)* gradDphi(1, 1) + gradDphi(1, 1) * gradDphi(0, 1), gradDphi(0, 2)* gradDphi(1, 1) + gradDphi(1, 2) * gradDphi(0, 1),
			gradDphi(0, 0)* gradDphi(1, 2) + gradDphi(1, 0) * gradDphi(0, 2), gradDphi(0, 1)* gradDphi(1, 2) + gradDphi(1, 1) * gradDphi(0, 2), gradDphi(0, 2)* gradDphi(1, 2) + gradDphi(1, 2) * gradDphi(0, 2);

		hessian->at(3) <<
			2 * gradDphi(1, 0) * gradDphi(1, 0), 2 * gradDphi(1, 1) * gradDphi(1, 0), 2 * gradDphi(1, 2) * gradDphi(1, 0),
			2 * gradDphi(1, 0) * gradDphi(1, 1), 2 * gradDphi(1, 1) * gradDphi(1, 1), 2 * gradDphi(1, 2) * gradDphi(1, 1),
			2 * gradDphi(1, 0) * gradDphi(1, 2), 2 * gradDphi(1, 1) * gradDphi(1, 2), 2 * gradDphi(1, 2) * gradDphi(1, 2);
	}

	return dphi * dphi.transpose();
}

std::vector<Eigen::Matrix2d> WTFShell::computeStretchingDensityFromQuad(int faceId, int quadId, std::vector<Eigen::MatrixXd >* deriv, std::vector<std::vector<Eigen::MatrixXd>>* hessian)
{
	std::vector<Eigen::Matrix2d> Mats;

	// compute some common terms
	Eigen::Vector2d da, dphi;
	Eigen::Vector3d gradA;
	Eigen::Matrix<double, 2, 3> gradDA, gradDphi;
	Eigen::Matrix<double, 3, 3> hessA;
	std::vector<Eigen::Matrix<double, 3, 3>> hessDA;

	double a = computeAmplitudesFromQuad(faceId, quadId, &da, deriv ? &gradA : NULL, deriv ? &gradDA : NULL, hessian ? &hessA : NULL, hessian ? &hessDA : NULL);
	// std::cout<<"amp : "<<a<<std::endl;
	dphi = computeDphi(faceId, deriv ? &gradDphi : NULL);
	Eigen::Matrix2d Ibar = _setup.abars[faceId];
	Eigen::Matrix2d I = firstFundamentalForm(_state.baseMesh, _state.basePos, faceId, NULL, NULL);
	Eigen::Matrix2d II = -_setup.sff->secondFundamentalForm(_state.baseMesh, _state.basePos, _state.baseEdgeDOFs, faceId, NULL, NULL);


	std::vector<Eigen::Matrix2d> grad_dadaTs;
	std::vector<Eigen::Matrix<double, 3, 3> > hess_dadaTs;
	Eigen::Matrix2d dadaTensor = computeDaDaTensor(faceId, quadId, (deriv || hessian) ? &grad_dadaTs : NULL, hessian ? &hess_dadaTs : NULL);

	std::vector<Eigen::Matrix2d> grad_dphidphiTs;
	std::vector<Eigen::Matrix<double, 3, 3> > hess_dphidphiTs;
	Eigen::Matrix2d dphidphiTensor = computeDphiDphiTensor(faceId, quadId, (deriv || hessian) ? &grad_dphidphiTs : NULL, hessian ? &hess_dphidphiTs : NULL);


	std::vector<Eigen::Matrix2d> grad_dphidaTs, grad_dadphiTs;
	std::vector<Eigen::Matrix<double, 6, 6>> hess_dphidaTs, hess_dadphiTs;

	Eigen::Matrix2d dadphiTensor = computeDaDphiTensor(faceId, quadId, (deriv || hessian) ? &grad_dadphiTs : NULL, hessian ? &hess_dadphiTs : NULL);
	Eigen::Matrix2d dphidaTensor = computeDphiDaTensor(faceId, quadId, (deriv || hessian) ? &grad_dphidaTs : NULL, hessian ? &hess_dphidaTs : NULL);


	std::vector<Eigen::Matrix2d> grad_dphidphiTSquares;
	std::vector<Eigen::Matrix<double, 3, 3> > hess_dphidphiTSquares;
	
	// compute the first term, namely the constant term in the Fourier series
	Eigen::Matrix2d M0 = Ibar.inverse() * (I + 0.5 * a * a * dphidphiTensor + 0.5 * dadaTensor - Ibar);

	Mats.push_back(M0);
	if (deriv)
	{
		Eigen::MatrixXd gradM0;
		gradM0.resize(4, 9);
		gradM0.setZero();
		Eigen::Matrix2d dM0dx;
		dM0dx.setZero();

		for (int i = 0; i < 3; i++)
		{
			dM0dx = Ibar.inverse() * (a * gradA(i) * dphidphiTensor + 0.5 * grad_dadaTs[i]);

			gradM0(0, i) = dM0dx(0, 0);
			gradM0(1, i) = dM0dx(0, 1);
			gradM0(2, i) = dM0dx(1, 0);
			gradM0(3, i) = dM0dx(1, 1);
		}

		for (int i = 0; i < 3; i++)
		{
			dM0dx = Ibar.inverse() * (0.5 * a * a * grad_dphidphiTs[i]);
			gradM0(0, 3 + i) = dM0dx(0, 0);
			gradM0(1, 3 + i) = dM0dx(0, 1);
			gradM0(2, 3 + i) = dM0dx(1, 0);
			gradM0(3, 3 + i) = dM0dx(1, 1);
		}

		
		deriv->push_back(gradM0);
	}

	if (hessian)
	{
		std::vector<Eigen::MatrixXd> hessM0;
		hessM0.resize(4);
		Eigen::Matrix2d IbarInv = Ibar.inverse();

		for (int i = 0; i < 4; i++)
		{
			hessM0[i].resize(9, 9);
			hessM0[i].setZero();
		}


		// hessian_aa part
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				hessM0[0](i, j) = 0.5 * IbarInv(0, 0) * ((2 * gradA(i) * gradA(j) + 2 * a * hessA(i, j)) * dphidphiTensor(0, 0) + hess_dadaTs[0](i, j));
				hessM0[0](i, j) += 0.5 * IbarInv(0, 1) * ((2 * gradA(i) * gradA(j) + 2 * a * hessA(i, j)) * dphidphiTensor(1, 0) + hess_dadaTs[2](i, j));

				hessM0[1](i, j) = 0.5 * IbarInv(0, 0) * ((2 * gradA(i) * gradA(j) + 2 * a * hessA(i, j)) * dphidphiTensor(0, 1) + hess_dadaTs[1](i, j));
				hessM0[1](i, j) += 0.5 * IbarInv(0, 1) * ((2 * gradA(i) * gradA(j) + 2 * a * hessA(i, j)) * dphidphiTensor(1, 1) + hess_dadaTs[3](i, j));

				hessM0[2](i, j) = 0.5 * IbarInv(1, 0) * ((2 * gradA(i) * gradA(j) + 2 * a * hessA(i, j)) * dphidphiTensor(0, 0) + hess_dadaTs[0](i, j));
				hessM0[2](i, j) += 0.5 * IbarInv(1, 1) * ((2 * gradA(i) * gradA(j) + 2 * a * hessA(i, j)) * dphidphiTensor(1, 0) + hess_dadaTs[1](i, j));

				hessM0[3](i, j) = 0.5 * IbarInv(1, 0) * ((2 * gradA(i) * gradA(j) + 2 * a * hessA(i, j)) * dphidphiTensor(0, 1) + hess_dadaTs[1](i, j));
				hessM0[3](i, j) += 0.5 * IbarInv(1, 1) * ((2 * gradA(i) * gradA(j) + 2 * a * hessA(i, j)) * dphidphiTensor(1, 1) + hess_dadaTs[3](i, j));
			}
		}

		// hessian_aphi part
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				hessM0[0](i, j + 3) = 0.5 * IbarInv(0, 0) * 2 * a * gradA(i) * grad_dphidphiTs[j](0, 0) + 0.5 * IbarInv(0, 1) * 2 * a * gradA(i) * grad_dphidphiTs[j](1, 0);
				hessM0[1](i, j + 3) = 0.5 * IbarInv(0, 0) * 2 * a * gradA(i) * grad_dphidphiTs[j](0, 1) + 0.5 * IbarInv(0, 1) * 2 * a * gradA(i) * grad_dphidphiTs[j](1, 1);
				hessM0[2](i, j + 3) = 0.5 * IbarInv(1, 0) * 2 * a * gradA(i) * grad_dphidphiTs[j](0, 0) + 0.5 * IbarInv(1, 1) * 2 * a * gradA(i) * grad_dphidphiTs[j](1, 0);
				hessM0[3](i, j + 3) = 0.5 * IbarInv(1, 0) * 2 * a * gradA(i) * grad_dphidphiTs[j](0, 1) + 0.5 * IbarInv(1, 1) * 2 * a * gradA(i) * grad_dphidphiTs[j](1, 1);
			}
		}

		// hessian_phia part
		for (int i = 0; i < 4; i++)
		{
			hessM0[i].block(3, 0, 3, 3) = hessM0[i].block(0, 3, 3, 3).transpose();
		}


		// hessian_phiphi part
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				hessM0[0](i + 3, j + 3) = 0.5 * IbarInv(0, 0) * a * a * hess_dphidphiTs[0](i, j) + 0.5 * IbarInv(0, 1) * a * a * hess_dphidphiTs[2](i, j);
				hessM0[1](i + 3, j + 3) = 0.5 * IbarInv(0, 0) * a * a * hess_dphidphiTs[1](i, j) + 0.5 * IbarInv(0, 1) * a * a * hess_dphidphiTs[3](i, j);
				hessM0[2](i + 3, j + 3) = 0.5 * IbarInv(1, 0) * a * a * hess_dphidphiTs[0](i, j) + 0.5 * IbarInv(1, 1) * a * a * hess_dphidphiTs[2](i, j);
				hessM0[3](i + 3, j + 3) = 0.5 * IbarInv(1, 0) * a * a * hess_dphidphiTs[1](i, j) + 0.5 * IbarInv(1, 1) * a * a * hess_dphidphiTs[3](i, j);

			}
		}

		hessian->push_back(hessM0);
	}

	// second term
	// formula without v1 in-plane correction tensileProj = -IbarInv * II
	Eigen::Matrix2d tensileProj;
	tensileProj.setZero();

	Eigen::Matrix2d A = Ibar.inverse() * (I - Ibar);
	Eigen::Vector2d evals;
	Eigen::Matrix2d evecs;

	Eigen::EigenSolver<Eigen::Matrix2d> solver;
	solver.compute(A);
	evals = solver.eigenvalues().real();
	evecs = solver.eigenvectors().real();


	double lambda1 = evals(0);
	double lambda2 = evals(1);

	Eigen::Vector2d compressDir;
	Eigen::Vector2d tensionDir;

	if (lambda1 <= 0 && lambda2 >= 0)
	{
		compressDir = evecs.col(0);
		tensionDir = evecs.col(1);
	}
	else if (lambda2 <= 0 && lambda1 >= 0)
	{
		compressDir = evecs.col(1);
		tensionDir = evecs.col(0);
	}
	else
	{
		tensionDir.setZero();
		tensileProj.setZero();
	}

	if(tensionDir.norm() > 0)
	{
		tensionDir = tensionDir / std::sqrt(tensionDir.dot(Ibar * tensionDir));
		double lameAlpha = _setup.YoungsModulus * _setup.PoissonsRatio / (1.0 - _setup.PoissonsRatio * _setup.PoissonsRatio);
		double lameBeta = _setup.YoungsModulus / 2.0 / (1.0 + _setup.PoissonsRatio);

		// tensileProj = std::sqrt(lameBeta * (lameAlpha + lameBeta)) / (lameAlpha / 2.0 + lameBeta) * II * tensionDir * tensionDir.transpose();
		Eigen::Matrix2d idMat;
		idMat.setIdentity();

		Eigen::Vector3d ru, rv;
		ru = _state.basePos.row(_state.baseMesh.faceVertex(faceId, 1)) - _state.basePos.row(_state.baseMesh.faceVertex(faceId, 0));
		rv = _state.basePos.row(_state.baseMesh.faceVertex(faceId, 2)) - _state.basePos.row(_state.baseMesh.faceVertex(faceId, 0));
		Eigen::Vector3d n = ru.cross(rv);
		n.normalize();
		Eigen::Vector3d tensileDir = tensionDir(0) * ru + tensionDir(1) * rv;
		ru.normalize();
		rv.normalize();
		tensileDir.normalize();
		double curvature = 0;
		for(int i = 0; i < 3; i++)
		{
			int vid = _state.baseMesh.faceVertex(faceId, i);
			Eigen::Vector3d pd1 = _state.PD1.row(vid);
			pd1 = pd1 - pd1.dot(n)*n;   // projection to the tangent plane
			
			double cosTheta = pd1.dot(tensileDir);
			
			curvature += _state.PV1(vid) * cosTheta * cosTheta + _state.PV2(vid) * (1 - cosTheta * cosTheta);
			
		}
		
		curvature = curvature / 3.0;

		tensileProj = std::sqrt(lameBeta / (lameAlpha + 2.0 * lameBeta)) * curvature * tensionDir.dot(I * tensionDir) * idMat;
	}
	else
	{
		tensileProj.setZero();
	}
	

	Eigen::Matrix2d M1 = 2 * a * tensileProj;

	Mats.push_back(M1);

	if (deriv)
	{
		Eigen::MatrixXd gradM1;
		gradM1.resize(4, 9);
		gradM1.setZero();
		Eigen::Matrix2d dM1dx;
		dM1dx.setZero();

		for (int i = 0; i < 3; i++)
		{
			dM1dx = (2.0 * tensileProj * gradA(i));

			gradM1(0, i) += dM1dx(0, 0);
			gradM1(1, i) += dM1dx(0, 1);
			gradM1(2, i) += dM1dx(1, 0);
			gradM1(3, i) += dM1dx(1, 1);
		}

		deriv->push_back(gradM1);
	}

	if (hessian)
	{
		std::vector<Eigen::MatrixXd> hessM1;
		hessM1.resize(4);
		for (int i = 0; i < 4; i++)
		{
			hessM1[i].resize(9, 9);
			hessM1[i].setZero();
		}
		   

		for(int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				Eigen::Matrix2d hessMat = (2.0 * hessA(i, j) * tensileProj);

				hessM1[0](i, j) += hessMat(0, 0);
				hessM1[1](i, j) += hessMat(0, 1);
				hessM1[2](i, j) += hessMat(1, 0);
				hessM1[3](i, j) += hessMat(1, 1);
			}

		hessian->push_back(hessM1);
	}

	// third term
	Eigen::Matrix2d M2 = 0.5 * Ibar.inverse() * dadaTensor;
	Mats.push_back(M2);

	if (deriv)
	{
		Eigen::MatrixXd gradM2;
		gradM2.resize(4, 9);
		gradM2.setZero();
		Eigen::Matrix2d dM2dx;
		dM2dx.setZero();

		for (int i = 0; i < 3; i++)
		{
			dM2dx = 0.5 * Ibar.inverse() * grad_dadaTs[i];

			gradM2(0, i) += dM2dx(0, 0);
			gradM2(1, i) += dM2dx(0, 1);
			gradM2(2, i) += dM2dx(1, 0);
			gradM2(3, i) += dM2dx(1, 1);
		}

		deriv->push_back(gradM2);
	}

	if (hessian)
	{
		std::vector<Eigen::MatrixXd > hessM2;
		hessM2.resize(4);
		Eigen::Matrix2d IbarInv = Ibar.inverse();

		for (int i = 0; i < 4; i++)
		{
			hessM2[i].resize(9, 9);
			hessM2[i].setZero();
		}


		// hessian_a part
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				hessM2[0](i, j) = 0.5 * IbarInv(0, 0) * hess_dadaTs[0](i, j);
				hessM2[0](i, j) += 0.5 * IbarInv(0, 1) * hess_dadaTs[2](i, j);

				hessM2[1](i, j) = 0.5 * IbarInv(0, 0) * hess_dadaTs[1](i, j);
				hessM2[1](i, j) += 0.5 * IbarInv(0, 1) * hess_dadaTs[3](i, j);

				hessM2[2](i, j) = 0.5 * IbarInv(1, 0) * hess_dadaTs[0](i, j);
				hessM2[2](i, j) += 0.5 * IbarInv(1, 1) * hess_dadaTs[1](i, j);

				hessM2[3](i, j) = 0.5 * IbarInv(1, 0) * hess_dadaTs[1](i, j);
				hessM2[3](i, j) += 0.5 * IbarInv(1, 1) * hess_dadaTs[3](i, j);

			}
		}

		hessian->push_back(hessM2);
	}

	// fourth one 
	Eigen::Matrix<double, 4, 6> d2phiDeriv;
	std::vector<Eigen::Matrix<double, 6, 6>> d2phiHessian;
	// Eigen::Matrix2d d2Phi = computeD2PhiPerface(faceId, (deriv || hessian) ? &d2phiDeriv : NULL, hessian ? &d2phiHessian : NULL);
	Eigen::Matrix2d M3;

	//if (_isUseD2phi)
	//    M3 = Ibar.inverse() * (-0.5 * a * (dphidaTensor + dadphiTensor) / 2.0 + 0.25 * a * a * d2Phi);
	//else
		M3 = Ibar.inverse() * (-0.5 * a * (dphidaTensor + dadphiTensor) / 2.0);

	Mats.push_back(M3);

	if (deriv)
	{
		Eigen::MatrixXd gradM3;
		gradM3.resize(4, 9);
		gradM3.setZero();
		Eigen::Matrix2d dM3dx;
		dM3dx.setZero();

		// a * dphi * da part
		for (int i = 0; i < 3; i++)
		{
			dM3dx = Ibar.inverse() * (-0.5 * gradA(i) * (dphidaTensor + dadphiTensor) / 2.0);
			gradM3(0, i) += dM3dx(0, 0);
			gradM3(1, i) += dM3dx(0, 1);
			gradM3(2, i) += dM3dx(1, 0);
			gradM3(3, i) += dM3dx(1, 1);
		}

		for (int i = 0; i < 6; i++)
		{
			Eigen::Matrix2d d2phiDerivMat;
			d2phiDerivMat << d2phiDeriv(0, i), d2phiDeriv(1, i), d2phiDeriv(2, i), d2phiDeriv(3, i);
			dM3dx = Ibar.inverse() * (-0.5 * a * (grad_dphidaTs[i] + grad_dadphiTs[i]) / 2.0);

			gradM3(0, i) += dM3dx(0, 0);
			gradM3(1, i) += dM3dx(0, 1);
			gradM3(2, i) += dM3dx(1, 0);
			gradM3(3, i) += dM3dx(1, 1);
		}

		// a^2 d2phi part

		deriv->push_back(gradM3);
	}


	if (hessian)
	{
		std::vector<Eigen::MatrixXd > hessM3;
		hessM3.resize(4);
		Eigen::Matrix2d IbarInv = Ibar.inverse();

		for (int i = 0; i < 4; i++)
		{
			hessM3[i].resize(9, 9);
			hessM3[i].setZero();
		}

		// a * dphi * da part

		// hessian_aa part
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				hessM3[0](i, j) += -0.5 * IbarInv(0, 0) / 2.0 * (gradA(i) * (grad_dphidaTs[j](0, 0) + grad_dadphiTs[j](0, 0)) + hessA(i, j) * (dphidaTensor(0, 0) + dadphiTensor(0, 0)) + gradA(j) * (grad_dphidaTs[i](0, 0) + grad_dadphiTs[i](0, 0)) + a * (hess_dadphiTs[0](i, j) + hess_dphidaTs[0](i, j)));
				hessM3[0](i, j) += -0.5 * IbarInv(0, 1) / 2.0 * (gradA(i) * (grad_dphidaTs[j](1, 0) + grad_dadphiTs[j](1, 0)) + hessA(i, j) * (dphidaTensor(1, 0) + dadphiTensor(1, 0)) + gradA(j) * (grad_dphidaTs[i](1, 0) + grad_dadphiTs[i](1, 0)) + a * (hess_dadphiTs[2](i, j) + hess_dphidaTs[2](i, j)));

				hessM3[1](i, j) += -0.5 * IbarInv(0, 0) / 2.0 * (gradA(i) * (grad_dphidaTs[j](0, 1) + grad_dadphiTs[j](0, 1)) + hessA(i, j) * (dphidaTensor(0, 1) + dadphiTensor(0, 1)) + gradA(j) * (grad_dphidaTs[i](0, 1) + grad_dadphiTs[i](0, 1)) + a * (hess_dadphiTs[1](i, j) + hess_dphidaTs[1](i, j)));
				hessM3[1](i, j) += -0.5 * IbarInv(0, 1) / 2.0 * (gradA(i) * (grad_dphidaTs[j](1, 1) + grad_dadphiTs[j](1, 1)) + hessA(i, j) * (dphidaTensor(1, 1) + dadphiTensor(1, 1)) + gradA(j) * (grad_dphidaTs[i](1, 1) + grad_dadphiTs[i](1, 1)) + a * (hess_dadphiTs[3](i, j) + hess_dphidaTs[3](i, j)));

				hessM3[2](i, j) += -0.5 * IbarInv(1, 0) / 2.0 * (gradA(i) * (grad_dphidaTs[j](0, 0) + grad_dadphiTs[j](0, 0)) + hessA(i, j) * (dphidaTensor(0, 0) + dadphiTensor(0, 0)) + gradA(j) * (grad_dphidaTs[i](0, 0) + grad_dadphiTs[i](0, 0)) + a * (hess_dadphiTs[0](i, j) + hess_dphidaTs[0](i, j)));
				hessM3[2](i, j) += -0.5 * IbarInv(1, 1) / 2.0 * (gradA(i) * (grad_dphidaTs[j](1, 0) + grad_dadphiTs[j](1, 0)) + hessA(i, j) * (dphidaTensor(1, 0) + dadphiTensor(1, 0)) + gradA(j) * (grad_dphidaTs[i](1, 0) + grad_dadphiTs[i](1, 0)) + a * (hess_dadphiTs[2](i, j) + hess_dphidaTs[2](i, j)));

				hessM3[3](i, j) += -0.5 * IbarInv(1, 0) / 2.0 * (gradA(i) * (grad_dphidaTs[j](0, 1) + grad_dadphiTs[j](0, 1)) + hessA(i, j) * (dphidaTensor(0, 1) + dadphiTensor(0, 1)) + gradA(j) * (grad_dphidaTs[i](0, 1) + grad_dadphiTs[i](0, 1)) + a * (hess_dadphiTs[1](i, j) + hess_dphidaTs[1](i, j)));
				hessM3[3](i, j) += -0.5 * IbarInv(1, 1) / 2.0 * (gradA(i) * (grad_dphidaTs[j](1, 1) + grad_dadphiTs[j](1, 1)) + hessA(i, j) * (dphidaTensor(1, 1) + dadphiTensor(1, 1)) + gradA(j) * (grad_dphidaTs[i](1, 1) + grad_dadphiTs[i](1, 1)) + a * (hess_dadphiTs[3](i, j) + hess_dphidaTs[3](i, j)));

			}
		}

		// hessian_aphi part
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				hessM3[0](i, 3 + j) = -0.5 * IbarInv(0, 0) / 2.0 * (gradA(i) * (grad_dphidaTs[3 + j](0, 0) + grad_dadphiTs[3 + j](0, 0)) + a * (hess_dadphiTs[0](i, 3 + j) + hess_dphidaTs[0](i, 3 + j)));
				hessM3[0](i, 3 + j) += -0.5 * IbarInv(0, 1) / 2.0 * (gradA(i) * (grad_dphidaTs[3 + j](1, 0) + grad_dadphiTs[3 + j](1, 0)) + a * (hess_dadphiTs[2](i, 3 + j) + hess_dphidaTs[2](i, 3 + j)));

				hessM3[1](i, 3 + j) = -0.5 * IbarInv(0, 0) / 2.0 * (gradA(i) * (grad_dphidaTs[3 + j](0, 1) + grad_dadphiTs[3 + j](0, 1)) + a * (hess_dadphiTs[1](i, 3 + j) + hess_dphidaTs[1](i, 3 + j)));
				hessM3[1](i, 3 + j) += -0.5 * IbarInv(0, 1) / 2.0 * (gradA(i) * (grad_dphidaTs[3 + j](1, 1) + grad_dadphiTs[3 + j](1, 1)) + a * (hess_dadphiTs[3](i, 3 + j) + hess_dphidaTs[3](i, 3 + j)));

				hessM3[2](i, 3 + j) = -0.5 * IbarInv(1, 0) / 2.0 * (gradA(i) * (grad_dphidaTs[3 + j](0, 0) + grad_dadphiTs[3 + j](0, 0)) + a * (hess_dadphiTs[0](i, 3 + j) + hess_dphidaTs[0](i, 3 + j)));
				hessM3[2](i, 3 + j) += -0.5 * IbarInv(1, 1) / 2.0 * (gradA(i) * (grad_dphidaTs[3 + j](1, 0) + grad_dadphiTs[3 + j](1, 0)) + a * (hess_dadphiTs[2](i, 3 + j) + hess_dphidaTs[2](i, 3 + j)));

				hessM3[3](i, 3 + j) = -0.5 * IbarInv(1, 0) / 2.0 * (gradA(i) * (grad_dphidaTs[3 + j](0, 1) + grad_dadphiTs[3 + j](0, 1)) + a * (hess_dadphiTs[1](i, 3 + j) + hess_dphidaTs[1](i, 3 + j)));
				hessM3[3](i, 3 + j) += -0.5 * IbarInv(1, 1) / 2.0 * (gradA(i) * (grad_dphidaTs[3 + j](1, 1) + grad_dadphiTs[3 + j](1, 1)) + a * (hess_dadphiTs[3](i, 3 + j) + hess_dphidaTs[3](i, 3 + j)));
			}
		}

		// hessian_phia part
		for (int i = 0; i < 4; i++)
		{
			hessM3[i].block(3, 0, 3, 3) = hessM3[i].block(0, 3, 3, 3).transpose();
		}


		// hessian_phiphi part
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				hessM3[0](3 + i, 3 + j) = -0.5 * IbarInv(0, 0) / 2.0 * (a * (hess_dadphiTs[0](3 + i, 3 + j) + hess_dphidaTs[0](3 + i, 3 + j)));
				hessM3[0](3 + i, 3 + j) += -0.5 * IbarInv(0, 1) / 2.0 * (a * (hess_dadphiTs[2](3 + i, 3 + j) + hess_dphidaTs[2](3 + i, 3 + j)));

				hessM3[1](3 + i, 3 + j) = -0.5 * IbarInv(0, 0) / 2.0 * (a * (hess_dadphiTs[1](3 + i, 3 + j) + hess_dphidaTs[1](3 + i, 3 + j)));
				hessM3[1](3 + i, 3 + j) += -0.5 * IbarInv(0, 1) / 2.0 * (a * (hess_dadphiTs[3](3 + i, 3 + j) + hess_dphidaTs[3](3 + i, 3 + j)));

				hessM3[2](3 + i, 3 + j) = -0.5 * IbarInv(1, 0) / 2.0 * (a * (hess_dadphiTs[0](3 + i, 3 + j) + hess_dphidaTs[0](3 + i, 3 + j)));
				hessM3[2](3 + i, 3 + j) += -0.5 * IbarInv(1, 1) / 2.0 * (a * (hess_dadphiTs[2](3 + i, 3 + j) + hess_dphidaTs[2](3 + i, 3 + j)));

				hessM3[3](3 + i, 3 + j) = -0.5 * IbarInv(1, 0) / 2.0 * (a * (hess_dadphiTs[1](3 + i, 3 + j) + hess_dphidaTs[1](3 + i, 3 + j)));
				hessM3[3](3 + i, 3 + j) += -0.5 * IbarInv(1, 1) / 2.0 * (a * (hess_dadphiTs[3](3 + i, 3 + j) + hess_dphidaTs[3](3 + i, 3 + j)));
			}
		}


		// a^2 * d2phi part


		hessian->push_back(hessM3);
	}
	return Mats;

}

std::vector<Eigen::Matrix2d> WTFShell::computeBendingDensityFromQuad(int faceId, int quadId, std::vector<Eigen::MatrixXd >* deriv, std::vector<std::vector<Eigen::MatrixXd>>* hessian)
{
	std::vector<Eigen::Matrix2d> Mats;

	// computer some common terms
	Eigen::Vector2d da, dphi;
	Eigen::Matrix2d d2a;
	Eigen::Vector3d gradA;
	Eigen::Matrix<double, 2, 3> gradDA, gradDphi;
	Eigen::Matrix<double, 4, 3> gradD2A;
	Eigen::Matrix<double, 3, 3> hessA;
	std::vector<Eigen::Matrix<double, 3, 3>> hessDA;

	double a = computeAmplitudesFromQuad(faceId, quadId, &da, (deriv || hessian) ? &gradA : NULL, deriv ? &gradDA : NULL, hessian ? &hessA : NULL, hessian ? &hessDA : NULL);
	dphi = computeDphi(faceId, deriv ? &gradDphi : NULL);

	Eigen::Matrix2d I = firstFundamentalForm(_state.baseMesh, _state.basePos, faceId, NULL, NULL);
	Eigen::Matrix2d II = -_setup.sff->secondFundamentalForm(_state.baseMesh, _state.basePos, _state.baseEdgeDOFs, faceId, NULL, NULL);
	Eigen::Matrix2d Ibar = _setup.abars[faceId];
	Eigen::Matrix2d IIbar = _setup.bbars[faceId];

	std::vector<Eigen::Matrix2d> grad_dphidphiTs;
	std::vector<Eigen::Matrix<double, 3, 3> > hess_dphidphiTs;
	Eigen::Matrix2d dphidphiTensor = computeDphiDphiTensor(faceId, quadId, (deriv || hessian) ? &grad_dphidphiTs : NULL, hessian ? &hess_dphidphiTs : NULL);
	
	// constant term
	// Eigen::Matrix2d M0 = Ibar.inverse() * (coeff * II - IIbar);
	Eigen::Matrix2d M0 = Ibar.inverse() * (II - IIbar);
	Mats.push_back(M0);

	if (deriv)
	{
		Eigen::MatrixXd gradM0;
		gradM0.resize(4, 6);
		gradM0.setZero();
		deriv->push_back(gradM0);
	}

	if (hessian)
	{
		std::vector<Eigen::MatrixXd > hessM0;
		hessM0.resize(4);
		for (int i = 0; i < 4; i++)
		{
			hessM0[i].resize(6, 6);
			hessM0[i].setZero();
		}

		hessian->push_back(hessM0);
	}



	// cos(phi) term

	// Eigen::Matrix2d M1 = Ibar.inverse() * coeff * (-a * (1 + 3.0 / 8.0 * adphiSquared) * dphidphiTensor);
	Eigen::Matrix2d M1 = Ibar.inverse() * (-a * dphidphiTensor);
	Mats.push_back(M1);

	if (deriv)
	{

		Eigen::MatrixXd gradM1;
		gradM1.resize(4, 6);
		gradM1.setZero();
		Eigen::Matrix2d dM1dx;
		dM1dx.setZero();

		for (int i = 0; i < 3; i++)
		{
			// dM1dx = Ibar.inverse() * coeff * (-gradA(i) * (1 + 3.0 / 8.0 * adphiSquared) * dphidphiTensor - a * 3.0 / 8.0 * grad_adphiSquared(i) * dphidphiTensor);
			dM1dx = Ibar.inverse() * (-gradA(i) * dphidphiTensor);

			gradM1(0, i) = dM1dx(0, 0);
			gradM1(1, i) = dM1dx(0, 1);
			gradM1(2, i) = dM1dx(1, 0);
			gradM1(3, i) = dM1dx(1, 1);
		}

		for (int i = 0; i < 3; i++)
		{
			// dM1dx = Ibar.inverse() * coeff * (-a * (1 + 3.0 / 8.0 * adphiSquared) * grad_dphidphiTs[i] - a * 3.0 / 8.0 * grad_adphiSquared(i + 3) * dphidphiTensor);
			dM1dx = Ibar.inverse() * (-a * grad_dphidphiTs[i]);

			gradM1(0, i + 3) = dM1dx(0, 0);
			gradM1(1, i + 3) = dM1dx(0, 1);
			gradM1(2, i + 3) = dM1dx(1, 0);
			gradM1(3, i + 3) = dM1dx(1, 1);
		}

		deriv->push_back(gradM1);

	}

	if (hessian)
	{
		std::vector<Eigen::MatrixXd > hessM1;
		hessM1.resize(4);
		Eigen::Matrix2d IbarInv = Ibar.inverse();

		Eigen::VectorXd fullGradA(6);
		fullGradA.setZero();
		fullGradA.segment(0, 3) = gradA;

		Eigen::Matrix<double, 6, 6> fullHessA;
		fullHessA.setZero();
		fullHessA.block(0, 0, 3, 3) = hessA;

		std::vector<Eigen::Matrix2d> fullGrad_dphidphiTs;
		for (int i = 0; i < 6; i++)
		{
			if (i < 3)
				fullGrad_dphidphiTs.push_back(Eigen::Matrix2d::Zero());
			else
				fullGrad_dphidphiTs.push_back(grad_dphidphiTs[i - 3]);
		}

		std::vector<Eigen::Matrix<double, 6, 6>> fullHess_dphidphiTs;
		for (int i = 0; i < hess_dphidphiTs.size(); i++)
		{
			Eigen::Matrix<double, 6, 6> mat = Eigen::Matrix<double, 6, 6>::Zero();
			mat.block(3, 3, 3, 3) = hess_dphidphiTs[i];
			fullHess_dphidphiTs.push_back(mat);
		}

		for (int i = 0; i < 4; i++)
		{
			hessM1[i].resize(6, 6);
			hessM1[i].setZero();
		}

		// hessian matrix part
		for (int i = 0; i < 6; i++)
			for (int j = 0; j < 6; j++)
			{

				Eigen::Matrix2d hessMat = (-fullHessA(i, j)) * dphidphiTensor;
				hessMat += (-fullGradA(i) ) * fullGrad_dphidphiTs[j];
				hessMat += (-fullGradA(j) ) * fullGrad_dphidphiTs[i];
				Eigen::Matrix2d hessdphidphiT;
				hessdphidphiT << fullHess_dphidphiTs[0](i, j), fullHess_dphidphiTs[1](i, j), fullHess_dphidphiTs[2](i, j), fullHess_dphidphiTs[3](i, j);
				hessMat += -a * hessdphidphiT;

				hessMat = IbarInv * hessMat;

				
				hessM1[0](i, j) += hessMat(0, 0);
				hessM1[1](i, j) += hessMat(0, 1);
				hessM1[2](i, j) += hessMat(1, 0);
				hessM1[3](i, j) += hessMat(1, 1);
			}



		hessian->push_back(hessM1);
	}

	return Mats;
}

double WTFShell::stretchingEnergyPerface(int faceId, Eigen::VectorXd* deriv, Eigen::MatrixXd* hessian)
{
	double lameAlpha = _setup.YoungsModulus * _setup.PoissonsRatio / (1.0 - _setup.PoissonsRatio * _setup.PoissonsRatio);
	double lameBeta = _setup.YoungsModulus / 2.0 / (1.0 + _setup.PoissonsRatio);
//    lameBeta = 0;

	Eigen::Matrix2d Ibar = _setup.abars[faceId];

	double energy = 0;
	double area = 0.5 * sqrt(Ibar.determinant());

	if (deriv)
	{
		deriv->resize(9);
		deriv->setZero();
	}

	if (hessian)
	{
		hessian->resize(9, 9);
		hessian->setZero();
	}

	for (int i = 0; i < _setup.quadPoints.size(); i++)
	{
		std::vector<Eigen::MatrixXd> matderivs;
		std::vector<std::vector<Eigen::MatrixXd > > mathessian;

		std::vector<Eigen::Matrix2d> Mats = computeStretchingDensityFromQuad(faceId, i, (deriv || hessian) ? &matderivs : NULL, hessian ? &mathessian : NULL);
		Eigen::VectorXd dphidaDeriv;

		for (int k = 0; k < Mats.size(); k++)
		{
			auto M = Mats[k];
			double factor = k == 0 ? 1 : 0.5;
			energy += factor * ((M * M).trace() * lameBeta + M.trace() * M.trace() * lameAlpha / 2.0) * _setup.quadPoints[i].weight * area * _setup.thickness / 4.0;
//            energy += factor * (M.trace() * M.trace() * lameAlpha / 2.0) * _setup.quadPoints[i].weight * area * _setup.thickness / 4.0;

			if (deriv)
			{
				Eigen::Matrix2d dM;
				for (int j = 0; j < 9; j++)
				{
					dM << matderivs[k](0, j), matderivs[k](1, j),
						matderivs[k](2, j), matderivs[k](3, j);
					deriv->coeffRef(j) += factor * (2.0 * (dM * M).trace() * lameBeta + 2.0 * dM.trace() * M.trace() * lameAlpha / 2.0) * _setup.quadPoints[i].weight * area * _setup.thickness / 4.0;
				}
			}

			if (hessian)
			{
				for (int m = 0; m < 9; m++)
					for (int n = 0; n < 9; n++)
					{
						hessian->coeffRef(m, n) += _setup.quadPoints[i].weight * area * _setup.thickness / 4.0 * factor * (lameBeta * (2 * matderivs[k](0, m) * matderivs[k](0, n) + 2 * Mats[k](0, 0) * mathessian[k][0](m, n)
							+ 2 * matderivs[k](3, m) * matderivs[k](3, n) + 2 * Mats[k](1, 1) * mathessian[k][3](m, n)));
						hessian->coeffRef(m, n) += _setup.quadPoints[i].weight * area * _setup.thickness / 4.0 * factor * (lameBeta * (2 * matderivs[k](1, m) * matderivs[k](2, n) + 2 * Mats[k](1, 0) * mathessian[k][1](m, n)
							+ 2 * matderivs[k](2, m) * matderivs[k](1, n) + 2 * Mats[k](0, 1) * mathessian[k][2](m, n)));
						hessian->coeffRef(m, n) += _setup.quadPoints[i].weight * area * _setup.thickness / 4.0 * factor * (0.5 * lameAlpha * (2 * (matderivs[k](0, m) + matderivs[k](3, m)) * (matderivs[k](0, n) + matderivs[k](3, n))
							+ 2 * (Mats[k](0, 0) + Mats[k](1, 1)) * (mathessian[k][0](m, n) + mathessian[k][3](m, n))));
					}

			}
		}
	}

	if (_isPosHess && hessian)
	{
		Eigen::MatrixXd posHess = lowRankApprox(*hessian);
		*hessian = posHess;
	}

	return energy;
}

double WTFShell::bendingEnergyPerface(int faceId, Eigen::VectorXd* deriv, Eigen::MatrixXd* hessian)
{
	double lameAlpha = _setup.YoungsModulus * _setup.PoissonsRatio / (1.0 - _setup.PoissonsRatio * _setup.PoissonsRatio);
	double lameBeta = _setup.YoungsModulus / 2.0 / (1.0 + _setup.PoissonsRatio);
	// lameBeta = 0;

	Eigen::Matrix2d Ibar = _setup.abars[faceId];

	double energy = 0;
	double area = 0.5 * sqrt(Ibar.determinant());

	if (deriv)
	{
		deriv->resize(6);
		deriv->setZero();
	}

	if (hessian)
	{
		hessian->resize(6, 6);
		hessian->setZero();
	}

	for (int i = 0; i < _setup.quadPoints.size(); i++)
	{
		std::vector<Eigen::MatrixXd > matderivs;
		std::vector<std::vector<Eigen::MatrixXd>> mathessian;

		std::vector<Eigen::Matrix2d> Mats = computeBendingDensityFromQuad(faceId, i, (deriv || hessian) ? &matderivs : NULL, hessian ? &mathessian : NULL);

		for (int k = 0; k < Mats.size(); k++)
		{
			auto M = Mats[k];
			double factor = k == 0 ? 1 : 0.5;
			energy += factor * ((M * M).trace() * lameBeta + M.trace() * M.trace() * lameAlpha / 2.0) * _setup.quadPoints[i].weight * area * std::pow(_setup.thickness, 3.0) / 12.0;
			// energy += factor * (M * M.transpose()).trace() * _setup.quadPoints[i].weight * area * std::pow(_setup.thickness, 2.0) / 3.0;

			if (deriv)
			{
				Eigen::Matrix2d dM;
				for (int j = 0; j < 6; j++)
				{
					dM << matderivs[k](0, j), matderivs[k](1, j),
						matderivs[k](2, j), matderivs[k](3, j);

					deriv->coeffRef(j) += factor * (2.0 * (dM * M).trace() * lameBeta + 2.0 * dM.trace() * M.trace() * lameAlpha / 2.0) * _setup.quadPoints[i].weight * area * std::pow(_setup.thickness, 3.0) / 12.0;
					// deriv->coeffRef(j) += factor * 2.0 * (dM * M.transpose()).trace() * _setup.quadPoints[i].weight * area * std::pow(_setup.thickness, 2.0) / 3.0;
				}
			}

			if (hessian)
			{
				for (int m = 0; m < 6; m++)
					for (int n = 0; n < 6; n++)
					{
						hessian->coeffRef(m, n) += _setup.quadPoints[i].weight * area * std::pow(_setup.thickness, 3.0) / 12.0 * factor * (lameBeta * (2 * matderivs[k](0, m) * matderivs[k](0, n) + 2 * Mats[k](0, 0) * mathessian[k][0](m, n)
							+ 2 * matderivs[k](3, m) * matderivs[k](3, n) + 2 * Mats[k](1, 1) * mathessian[k][3](m, n)));
						hessian->coeffRef(m, n) += _setup.quadPoints[i].weight * area * std::pow(_setup.thickness, 3.0) / 12.0 * factor * (lameBeta * (2 * matderivs[k](1, m) * matderivs[k](2, n) + 2 * Mats[k](1, 0) * mathessian[k][1](m, n)
							+ 2 * matderivs[k](2, m) * matderivs[k](1, n) + 2 * Mats[k](0, 1) * mathessian[k][2](m, n)));
						hessian->coeffRef(m, n) += _setup.quadPoints[i].weight * area * std::pow(_setup.thickness, 3.0) / 12.0 * factor * (0.5 * lameAlpha * (2 * (matderivs[k](0, m) + matderivs[k](3, m)) * (matderivs[k](0, n) + matderivs[k](3, n))
							+ 2 * (Mats[k](0, 0) + Mats[k](1, 1)) * (mathessian[k][0](m, n) + mathessian[k][3](m, n))));
					}
			}
		}
	}

	if (_isPosHess && hessian)
	{
		Eigen::MatrixXd posHess = lowRankApprox(*hessian);
		*hessian = posHess;
	}

	return energy;
}

double WTFShell::stretchingEnergy(Eigen::VectorXd* deriv, Eigen::SparseMatrix<double>* hessian)
{
	int nfaces = _state.baseMesh.nFaces();
	int nverts = _state.basePos.rows();
	int nedges = _state.baseMesh.nEdges();
	double energy = 0;

	if (deriv)
	{
		deriv->resize(nverts + nedges);
		deriv->setZero();
	}

	std::vector<Eigen::Triplet<double> > hessianT;
	auto energies = std::vector<double>(nfaces);
	auto dStretchings = std::vector<Eigen::VectorXd>(nfaces);
	auto hStretchings = std::vector<Eigen::MatrixXd>(nfaces);

	if (_isParallel)
	{
		auto computeStretching = [&](const tbb::blocked_range<uint32_t>& range)
		{
			for (uint32_t i = range.begin(); i < range.end(); ++i)
				energies[i] += stretchingEnergyPerface(i, deriv ? &dStretchings[i] : NULL, hessian ? &hStretchings[i] : NULL);
		};

		tbb::blocked_range<uint32_t> rangex(0u, (uint32_t)nfaces, GRAIN_SIZE);
		tbb::parallel_for(rangex, computeStretching);
	}
	else
	{
		for (int i = 0; i < nfaces; i++)
			energies[i] += stretchingEnergyPerface(i, deriv ? &dStretchings[i] : NULL, hessian ? &hStretchings[i] : NULL);
	}

	

	for (int i = 0; i < nfaces; i++)
	{
		energy += energies[i];
		if (deriv)
		{
			for (int j = 0; j < 3; j++)
			{
				int vid = _state.baseMesh.faceVertex(i, j);
				int eid = _state.baseMesh.faceEdge(i, j);
				int oppid = _state.baseMesh.vertexOppositeFaceEdge(i, j);

				deriv->coeffRef(vid) += dStretchings[i](j);
				deriv->coeffRef(nverts + eid) += dStretchings[i](3 + j);
			}
		}

		if (hessian)
		{
			for (int j = 0; j < 3; j++)
			{
				int vid = _state.baseMesh.faceVertex(i, j);
				int eid = _state.baseMesh.faceEdge(i, j);
				int oppid = _state.baseMesh.vertexOppositeFaceEdge(i, j);

				for (int k = 0; k < 3; k++)
				{
					int vid1 = _state.baseMesh.faceVertex(i, k);
					int eid1 = _state.baseMesh.faceEdge(i, k);
				
					hessianT.push_back(Eigen::Triplet<double>(vid, vid1, hStretchings[i](j, k)));
					hessianT.push_back(Eigen::Triplet<double>(vid, eid1 + nverts, hStretchings[i](j, 3 + k)));

					hessianT.push_back(Eigen::Triplet<double>(eid + nverts, vid1, hStretchings[i](3 + j, k)));
					hessianT.push_back(Eigen::Triplet<double>(eid + nverts, eid1 + nverts, hStretchings[i](3 + j, 3 + k)));

				}

			}
		}
	}


	if (hessian)
	{
		hessian->resize(nverts + nedges, nverts + nedges);
		hessian->setFromTriplets(hessianT.begin(), hessianT.end());
	}

	return energy;
}

double WTFShell::bendingEnergy(Eigen::VectorXd* deriv, Eigen::SparseMatrix<double>* hessian)
{
	int nfaces = _state.baseMesh.nFaces();
	int nverts = _state.basePos.rows();
	int nedges = _state.baseMesh.nEdges();

	double energy = 0;
	if (deriv)
	{
		deriv->resize(nverts + nedges);
		deriv->setZero();
	}

	std::vector<Eigen::Triplet<double> > hessianT;

	auto energies = std::vector<double>(nfaces);
	auto dBendings = std::vector<Eigen::VectorXd>(nfaces);
	auto hBendings = std::vector<Eigen::MatrixXd>(nfaces);

	if (_isParallel)
	{
		auto computeBending = [&](const tbb::blocked_range<uint32_t>& range)
		{
			for (uint32_t i = range.begin(); i < range.end(); ++i)
				energies[i] += bendingEnergyPerface(i, deriv ? &dBendings[i] : NULL, hessian ? &hBendings[i] : NULL);
		};

		tbb::blocked_range<uint32_t> rangex(0u, (uint32_t)nfaces, GRAIN_SIZE);
		tbb::parallel_for(rangex, computeBending);
	}
	else
	{
		for (int i = 0; i < nfaces; i++)
			energies[i] += bendingEnergyPerface(i, deriv ? &dBendings[i] : NULL, hessian ? &hBendings[i] : NULL);
	}
	

	for (int i = 0; i < nfaces; i++)
	{
		energy += energies[i];
		if (deriv)
		{
			for (int j = 0; j < 3; j++)
			{
				int vid = _state.baseMesh.faceVertex(i, j);
				int eid = _state.baseMesh.faceEdge(i, j);

				deriv->coeffRef(vid) += dBendings[i](j);
				deriv->coeffRef(nverts + eid) += dBendings[i](3 + j);
			}
		}

		if (hessian)
		{
			for (int j = 0; j < 3; j++)
			{
				int vid = _state.baseMesh.faceVertex(i, j);
				int eid = _state.baseMesh.faceEdge(i, j);
				int oppid = _state.baseMesh.vertexOppositeFaceEdge(i, j);

				for (int k = 0; k < 3; k++)
				{
					int vid1 = _state.baseMesh.faceVertex(i, k);
					int eid1 = _state.baseMesh.faceEdge(i, k);

					int oppid1 = _state.baseMesh.vertexOppositeFaceEdge(i, k);

					hessianT.push_back(Eigen::Triplet<double>(vid, vid1, hBendings[i](j, k)));
					hessianT.push_back(Eigen::Triplet<double>(vid, eid1 + nverts, hBendings[i](j, 3 + k)));

					hessianT.push_back(Eigen::Triplet<double>(eid + nverts, vid1, hBendings[i](3 + j, k)));
					hessianT.push_back(Eigen::Triplet<double>(eid + nverts, eid1 + nverts, hBendings[i](3 + j, 3 + k)));

				}

			}
		}
	}

	if (hessian)
	{
		hessian->resize(nverts + nedges, nverts + nedges);
		hessian->setFromTriplets(hessianT.begin(), hessianT.end());
	}

	return energy;
}

double WTFShell::elasticReducedEnergy(Eigen::VectorXd* deriv, Eigen::SparseMatrix<double>* hessian)
{
	Eigen::VectorXd dEs, dEb, dEla, dEp;
	Eigen::SparseMatrix<double> hEs, hEb;
	double Es, Eb, Ela, Ep;
	Es = stretchingEnergy(deriv ? &dEs : NULL, hessian ? &hEs : NULL);
	Eb = bendingEnergy(deriv ? &dEb : NULL, hessian ? &hEb : NULL);

	if (deriv)
	{
		*deriv = dEs + dEb;
	}

	if (hessian)
	{
		*hessian = hEs + hEb;
	}

	return Es + Eb;
}

Eigen::VectorXd WTFShell::computeStretchingEnergyTermByTerm()
{
	double lameAlpha = _setup.YoungsModulus * _setup.PoissonsRatio / (1.0 - _setup.PoissonsRatio * _setup.PoissonsRatio);
	double lameBeta = _setup.YoungsModulus / 2.0 / (1.0 + _setup.PoissonsRatio);
	//    lameBeta = 0;

	Eigen::VectorXd stretchingEnergyTerms(5);
	stretchingEnergyTerms.setZero();

	for (int faceId = 0; faceId < _state.baseMesh.nFaces(); faceId++)
	{

		Eigen::Matrix2d Ibar = _setup.abars[faceId];

		double area = 0.5 * sqrt(Ibar.determinant());

		for (int i = 0; i < _setup.quadPoints.size(); i++)
		{
			std::vector<Eigen::MatrixXd> matderivs;
			std::vector<std::vector<Eigen::MatrixXd > > mathessian;

			std::vector<Eigen::Matrix2d> Mats = computeStretchingDensityFromQuad(faceId, i, NULL, NULL);
			Eigen::VectorXd dphidaDeriv;

			for (int k = 0; k < Mats.size(); k++)
			{
				double energy = 0;
				auto M = Mats[k];
				double factor = k == 0 ? 1 : 0.5;
				stretchingEnergyTerms(k) += factor * ((M * M).trace() * lameBeta + M.trace() * M.trace() * lameAlpha / 2.0) * _setup.quadPoints[i].weight * area * _setup.thickness / 4.0;
			}
		}
	}

	return stretchingEnergyTerms;
}

Eigen::VectorXd WTFShell::computeBendingEnergyTermByTerm()
{
	double lameAlpha = _setup.YoungsModulus * _setup.PoissonsRatio / (1.0 - _setup.PoissonsRatio * _setup.PoissonsRatio);
	double lameBeta = _setup.YoungsModulus / 2.0 / (1.0 + _setup.PoissonsRatio);
	//    lameBeta = 0;

	Eigen::VectorXd bendingEnergyTerms(4);
	bendingEnergyTerms.setZero();

	for (int faceId = 0; faceId < _state.baseMesh.nFaces(); faceId++)
	{

		Eigen::Matrix2d Ibar = _setup.abars[faceId];

		double area = 0.5 * sqrt(Ibar.determinant());

		for (int i = 0; i < _setup.quadPoints.size(); i++)
		{
			std::vector<Eigen::MatrixXd> matderivs;
			std::vector<std::vector<Eigen::MatrixXd > > mathessian;

			std::vector<Eigen::Matrix2d> Mats = computeBendingDensityFromQuad(faceId, i, NULL, NULL);
			Eigen::VectorXd dphidaDeriv;

			for (int k = 0; k < Mats.size(); k++)
			{
				double energy = 0;
				auto M = Mats[k];
				double factor = k == 0 ? 1 : 0.5;
				bendingEnergyTerms(k) += factor * ((M * M).trace() * lameBeta + M.trace() * M.trace() * lameAlpha / 2.0) * _setup.quadPoints[i].weight * area * std::pow(_setup.thickness, 3.0) / 12.0;
			}
		}
	}

	return bendingEnergyTerms;
}


//////////////////////////////////////// Test Functions ///////////////////////////////////////////////

void WTFShell::testAmplitudesFromQuad(int faceId, int quadId)
{
	auto backupSetup = _setup;
	auto backupState = _state;

	Eigen::Vector2d da;
	Eigen::Vector3d gradA;
	Eigen::Matrix<double, 2, 3> gradDA;

	Eigen::Matrix<double, 3, 3> hessA;
	std::vector<Eigen::Matrix<double, 3, 3>> hessDA;

	double a = computeAmplitudesFromQuad(faceId, quadId, &da, &gradA, &gradDA, &hessA, &hessDA);

	std::cout << "Test amplitude, face id: " << faceId << ", quad Id: " << quadId << ", u: " << _setup.quadPoints[quadId].u << ", v: " << _setup.quadPoints[quadId].v << std::endl;

	std::cout << std::endl << "test da: " << std::endl;
	Eigen::Vector2d dir = Eigen::Vector2d::Random();
	dir.normalize();
	for (int i = 3; i < 10; i++)
	{
		double eps = std::pow(10, -i);
		_setup.quadPoints[quadId].u = backupSetup.quadPoints[quadId].u + eps * dir(0);
		_setup.quadPoints[quadId].v = backupSetup.quadPoints[quadId].v + eps * dir(1);
		double a1 = computeAmplitudesFromQuad(faceId, quadId, NULL, NULL, NULL, NULL, NULL);
		std::cout << "eps: " << eps << std::endl;
		std::cout << "finite difference: " << (a1 - a) / eps << std::endl;
		std::cout << "directional derivative: " << da.dot(dir) << std::endl;
		std::cout << "error: " << std::abs((a1 - a) / eps - da.dot(dir)) << std::endl;
	}
	_setup = backupSetup;

	std::cout << std::endl << "test grad A: " << std::endl;
	Eigen::Vector3d dir1 = Eigen::Vector3d::Random();
	dir1.normalize();
	for (int i = 3; i < 10; i++)
	{
		double eps = std::pow(10, -i);
		for (int j = 0; j < 3; j++)
		{
			int vid = _state.baseMesh.faceVertex(faceId, j);
			_state.amplitude(vid) = backupState.amplitude(vid) + eps * dir1(j);
		}

		double a1 = computeAmplitudesFromQuad(faceId, quadId, NULL, NULL, NULL, NULL, NULL);
		std::cout << "eps: " << eps << std::endl;
		std::cout << "finite difference: " << (a1 - a) / eps << std::endl;
		std::cout << "directional derivative: " << gradA.dot(dir1) << std::endl;
		std::cout << "error: " << std::abs((a1 - a) / eps - gradA.dot(dir1)) << std::endl;
	}
	_state = backupState;

	std::cout << std::endl << "test grad dA: " << std::endl;
	for (int i = 3; i < 10; i++)
	{
		double eps = std::pow(10, -i);
		for (int j = 0; j < 3; j++)
		{
			int vid = _state.baseMesh.faceVertex(faceId, j);
			_state.amplitude(vid) = backupState.amplitude(vid) + eps * dir1(j);
		}

		Eigen::Vector2d da1;
		double a1 = computeAmplitudesFromQuad(faceId, quadId, &da1, NULL, NULL, NULL, NULL);
		std::cout << "eps: " << eps << std::endl;
		std::cout << "finite difference: " << ((da1 - da) / eps).transpose() << std::endl;
		std::cout << "directional derivative: " << (gradDA * dir1).transpose() << std::endl;
		std::cout << "error: " << ((da1 - da) / eps - gradDA * dir1).norm() << std::endl;
	}
	_state = backupState;



}

void WTFShell::testDphi(int faceId)
{
	Eigen::Matrix<double, 2, 3> gradDphi;
	Eigen::Vector2d dphi = computeDphi(faceId, &gradDphi);
	Eigen::Vector3d dir = Eigen::Vector3d::Random();
	dir.normalize();
	auto backupState = _state;
	std::cout << std::endl << "test grad dphi: " << std::endl;
	for (int i = 3; i < 10; i++)
	{
		double eps = std::pow(10, -i);
		for (int j = 0; j < 3; j++)
		{
			int eid = _state.baseMesh.faceEdge(faceId, j);
			_state.dphi(eid) = backupState.dphi(eid) + eps * dir(j);
		}

		Eigen::Vector2d dphi1 = computeDphi(faceId, NULL);
		std::cout << "eps: " << eps << std::endl;
		std::cout << "finite difference: " << ((dphi1 - dphi) / eps).transpose() << std::endl;
		std::cout << "directional derivative: " << (gradDphi * dir).transpose() << std::endl;
		std::cout << "error: " << ((dphi1 - dphi) / eps - gradDphi * dir).norm() << std::endl;
	}
	_state = backupState;
}

void WTFShell::testStretchingDensityFromQuad(int faceId, int quadId)
{
	std::vector<Eigen::MatrixXd > deriv;
	std::vector<std::vector<Eigen::MatrixXd>> hessian;
	std::vector<Eigen::Matrix2d> mats = computeStretchingDensityFromQuad(faceId, quadId, &deriv, &hessian);

	std::cout << "Test stretching density, face id: " << faceId << ", quad Id: " << quadId << ", u: " << _setup.quadPoints[quadId].u << ", v: " << _setup.quadPoints[quadId].v << std::endl;

	Eigen::Matrix<double, 9, 1> dir = Eigen::Matrix<double, 9, 1>::Random();
	dir.normalize();
	auto backupState = _state;

	for (int k = 0; k < mats.size(); k++)
	{
		std::cout << std::endl << k << "-th metric: " << std::endl;
		for (int i = 3; i < 10; i++)
		{
			double eps = std::pow(10, -i);
			for (int j = 0; j < 3; j++)
			{
				int vid = _state.baseMesh.faceVertex(faceId, j);
				int eid = _state.baseMesh.faceEdge(faceId, j);

				_state.amplitude(vid) = backupState.amplitude(vid) + eps * dir(j);
				_state.dphi(eid) = backupState.dphi(eid) + eps * dir(3 + j);
			}

			for (int j = 0; j < 3; j++)
			{
				/* int oppid = _state.baseMesh.vertexOppositeFaceEdge(faceId, j);
				 if (oppid != -1)
					 _state.phi(oppid) = backupState.phi(oppid) + eps * dir(6 + j);
				 else*/
				dir(6 + j) = 0;
			}

			std::vector<Eigen::MatrixXd > deriv1;
			std::vector<Eigen::Matrix2d> mats1 = computeStretchingDensityFromQuad(faceId, quadId, &deriv1, NULL);
			Eigen::Vector4d matVec, matVec1;
			matVec << mats[k](0, 0), mats[k](0, 1), mats[k](1, 0), mats[k](1, 1);
			matVec1 << mats1[k](0, 0), mats1[k](0, 1), mats1[k](1, 0), mats1[k](1, 1);

			std::cout << std::endl << "eps: " << eps << std::endl;
			std::cout << "energy - derivative: " << std::endl;
			std::cout << "finite differece: " << ((matVec1 - matVec) / eps).transpose() << std::endl;
			std::cout << "directional derivative: " << (deriv[k] * dir).transpose() << std::endl;
			std::cout << "error: " << ((matVec1 - matVec) / eps - deriv[k] * dir).norm() << std::endl;

			std::cout << std::endl << "derivative - hessian: " << std::endl;
			for (int j = 0; j < 4; j++)
			{
				std::cout << j << "-th entry: " << std::endl;
				std::cout << "finite diffference: " << (deriv1[k].row(j) - deriv[k].row(j)) / eps << std::endl;
				std::cout << "directional derivative: " << (hessian[k][j] * dir).transpose() << std::endl;
				std::cout << "error: " << ((hessian[k][j] * dir).transpose() - (deriv1[k].row(j) - deriv[k].row(j)) / eps).norm() << std::endl;
			}
		}
	}
	_state = backupState;
}

void WTFShell::testBendingDensityFromQuad(int faceId, int quadId)
{
	std::vector<Eigen::MatrixXd > deriv;
	std::vector<std::vector<Eigen::MatrixXd>> hessian;
	std::vector<Eigen::Matrix2d> mats = computeBendingDensityFromQuad(faceId, quadId, &deriv, &hessian);

	std::cout << "Test bending density, face id: " << faceId << ", quad Id: " << quadId << ", u: " << _setup.quadPoints[quadId].u << ", v: " << _setup.quadPoints[quadId].v << std::endl;

	Eigen::Matrix<double, 6, 1> dir = Eigen::Matrix<double, 6, 1>::Random();
	dir.normalize();
	auto backupState = _state;

	for (int k = 0; k < mats.size(); k++)
	{
		std::cout << std::endl << "test " << k << "-th term" << std::endl;
		for (int i = 3; i < 10; i++)
		{
			double eps = std::pow(10, -i);
			for (int j = 0; j < 3; j++)
			{
				int vid = _state.baseMesh.faceVertex(faceId, j);
				int eid = _state.baseMesh.faceEdge(faceId, j);

				_state.amplitude(vid) = backupState.amplitude(vid) + eps * dir(j);
				_state.dphi(eid) = backupState.dphi(eid) + eps * dir(3 + j);
			}

			std::vector<Eigen::MatrixXd > deriv1;
			std::vector<Eigen::Matrix2d> mats1 = computeBendingDensityFromQuad(faceId, quadId, &deriv1, NULL);

			Eigen::Vector4d matVec, matVec1;
			matVec << mats[k](0, 0), mats[k](0, 1), mats[k](1, 0), mats[k](1, 1);
			matVec1 << mats1[k](0, 0), mats1[k](0, 1), mats1[k](1, 0), mats1[k](1, 1);

			std::cout << std::endl << "eps: " << eps << std::endl;
			std::cout << "energy - derivative: " << std::endl;
			std::cout << "finite differece: " << ((matVec1 - matVec) / eps).transpose() << std::endl;
			std::cout << "directional derivative: " << (deriv[k] * dir).transpose() << std::endl;
			std::cout << "error: " << ((matVec1 - matVec) / eps - deriv[k] * dir).norm() << std::endl;

			std::cout << std::endl << "derivative - hessian: " << std::endl;
			for (int j = 0; j < 4; j++)
			{
				std::cout << j << "-th entry: " << std::endl;
				std::cout << "finite diffference: " << (deriv1[k].row(j) - deriv[k].row(j)) / eps << std::endl;
				std::cout << "directional derivative: " << (hessian[k][j] * dir).transpose() << std::endl;
				std::cout << "error: " << ((hessian[k][j] * dir).transpose() - (deriv1[k].row(j) - deriv[k].row(j)) / eps).norm() << std::endl;
			}
		}
	}
	_state = backupState;
}

void WTFShell::testStretchingEnergyPerface(int faceId)
{
	auto backupPosHess = _isPosHess;
	_isPosHess = false;

	Eigen::VectorXd deriv;
	Eigen::MatrixXd hessian;
	double energy = stretchingEnergyPerface(faceId, &deriv, &hessian);

	std::cout << std::endl << "Test stretching, face id: " << faceId << std::endl;

	Eigen::Matrix<double, 9, 1> dir = Eigen::Matrix<double, 9, 1>::Random();
	dir.normalize();
	auto backupState = _state;
	

	for (int i = 3; i < 10; i++)
	{
		double eps = std::pow(10, -i);
		for (int j = 0; j < 3; j++)
		{
			int vid = _state.baseMesh.faceVertex(faceId, j);
			int eid = _state.baseMesh.faceEdge(faceId, j);

			_state.amplitude(vid) = backupState.amplitude(vid) + eps * dir(j);
			_state.dphi(eid) = backupState.dphi(eid) + eps * dir(3 + j);
		}
		for (int j = 0; j < 3; j++)
		{
			/* int oppid = _state.baseMesh.vertexOppositeFaceEdge(faceId, j);
			 if (oppid != -1)
				 _state.phi(oppid) = backupState.phi(oppid) + eps * dir(6 + j);
			 else*/
			dir(6 + j) = 0;
		}

		Eigen::VectorXd deriv1;
		double energy1 = stretchingEnergyPerface(faceId, &deriv1, NULL);
		std::cout << std::endl << "eps: " << eps << std::endl;
		std::cout << std::endl << "energy - derivative: " << std::endl;
		std::cout << "finite difference: " << (energy1 - energy) / eps << std::endl;
		std::cout << "directional derivative: " << deriv.dot(dir) << std::endl;
		std::cout << "error: " << std::abs((energy1 - energy) / eps - deriv.dot(dir)) << std::endl;

		std::cout << std::endl << "derivative - hessian: " << std::endl;
		std::cout << "finite difference: " << (deriv1 - deriv).transpose() / eps << std::endl;
		std::cout << "directional derivative: " << (hessian * dir).transpose() << std::endl;
		std::cout << "error: " << ((deriv1 - deriv) / eps - hessian * dir).norm() << std::endl;
	}
	_state = backupState;
	_isPosHess = backupPosHess;

}

void WTFShell::testBendingEnergyPerface(int faceId)
{
	auto backupPosHess = _isPosHess;
	_isPosHess = false;

	Eigen::VectorXd deriv;
	Eigen::MatrixXd hessian;
	double energy = bendingEnergyPerface(faceId, &deriv, &hessian);

	std::cout << std::endl << "Test bending, face id: " << faceId << std::endl;

	Eigen::Matrix<double, 6, 1> dir = Eigen::Matrix<double, 6, 1>::Random();
	dir.normalize();
	auto backupState = _state;
	

	for (int i = 3; i < 10; i++)
	{
		double eps = std::pow(10, -i);
		for (int j = 0; j < 3; j++)
		{
			int vid = _state.baseMesh.faceVertex(faceId, j);
			int eid = _state.baseMesh.faceEdge(faceId, j);

			_state.amplitude(vid) = backupState.amplitude(vid) + eps * dir(j);
			_state.dphi(eid) = backupState.dphi(eid) + eps * dir(3 + j);
		}
		Eigen::VectorXd deriv1;
		double energy1 = bendingEnergyPerface(faceId, &deriv1, NULL);

		std::cout << std::endl << "eps: " << eps << std::endl;
		std::cout << std::endl << "energy - derivative: " << std::endl;
		std::cout << "finite difference: " << (energy1 - energy) / eps << std::endl;
		std::cout << "directional derivative: " << deriv.dot(dir) << std::endl;
		std::cout << "error: " << std::abs((energy1 - energy) / eps - deriv.dot(dir)) << std::endl;

		std::cout << std::endl << "derivative - hessian: " << std::endl;
		std::cout << "finite difference: " << (deriv1 - deriv).transpose() / eps << std::endl;
		std::cout << "directional derivative: " << (hessian * dir).transpose() << std::endl;
		std::cout << "error: " << ((deriv1 - deriv) / eps - hessian * dir).norm() << std::endl;
	}
	_state = backupState;
	_isPosHess = backupPosHess;
}

void WTFShell::testStretchingEnergy()
{
	auto backupPosHess = _isPosHess;
	_isPosHess = false;

	Eigen::VectorXd deriv;
	Eigen::SparseMatrix<double> hessian;

	double energy = stretchingEnergy(&deriv, &hessian);

	std::cout.precision(10);
	std::cout << std::endl << "Test stretching: " << std::endl;

	Eigen::VectorXd dir = Eigen::VectorXd::Random(deriv.size());
	dir.normalize();
	auto backupState = _state;

	int nverts = _state.basePos.rows();
	int nedges = _state.baseMesh.nEdges();

	std::cout << "verts: " << nverts << std::endl;
	std::cout << "edges: " << nedges << std::endl;
	std::cout << "energy: " << energy << std::endl;

	for (int i = 3; i < 10; i++)
	{
		double eps = std::pow(10, -i);
		for (int j = 0; j < nverts; j++)
		{
			_state.amplitude(j) = backupState.amplitude(j) + eps * dir(j);
		}
		for (int j = 0; j < nedges; j++)
		{
			_state.dphi(j) = backupState.dphi(j) + eps * dir(j + nverts);
		}

		Eigen::VectorXd deriv1;
		double energy1 = stretchingEnergy(&deriv1, NULL);
		//std::cout<<energy1<<std::endl;
		std::cout << std::endl << "eps: " << eps << std::endl;
		std::cout << std::endl << "energy - derivative: " << std::endl;
		std::cout << "finite difference: " << (energy1 - energy) / eps << std::endl;
		std::cout << "directional derivative: " << deriv.dot(dir) << std::endl;
		std::cout << "error: " << std::abs((energy1 - energy) / eps - deriv.dot(dir)) << std::endl;

		std::cout << std::endl << "derivative - hessian: " << std::endl;
		std::cout << "finite difference: " << (deriv1 - deriv).norm() / eps << std::endl;
		std::cout << "directional derivative: " << (hessian * dir).norm() << std::endl;
		std::cout << "error: " << ((deriv1 - deriv) / eps - hessian * dir).norm() << std::endl;
	}
	_state = backupState;
	_isPosHess = backupPosHess;

}

void WTFShell::testBendingEnergy()
{
	auto backupPosHess = _isPosHess;
	_isPosHess = false;

	Eigen::VectorXd deriv;
	Eigen::SparseMatrix<double> hessian;

	double energy = bendingEnergy(&deriv, &hessian);

	std::cout.precision(10);
	std::cout << std::endl << "Test bending: " << std::endl;

	Eigen::VectorXd dir = Eigen::VectorXd::Random(deriv.size());
	dir.normalize();
	auto backupState = _state;

	int nverts = _state.basePos.rows();
	int nedges = _state.baseMesh.nEdges();

	std::cout << "verts: " << nverts << std::endl;
	std::cout << "edges: " << nedges << std::endl;
	std::cout << "energy: " << energy << std::endl;

	for (int i = 3; i < 10; i++)
	{
		double eps = std::pow(10, -i);
		for (int j = 0; j < nverts; j++)
		{
			_state.amplitude(j) = backupState.amplitude(j) + eps * dir(j);
		}
		for (int j = 0; j < nedges; j++)
		{
			_state.dphi(j) = backupState.dphi(j) + eps * dir(j + nverts);
		}

		Eigen::VectorXd deriv1;
		double energy1 = bendingEnergy(&deriv1, NULL);
		//std::cout<<energy1<<std::endl;
		std::cout << std::endl << "eps: " << eps << std::endl;
		std::cout << std::endl << "energy - derivative: " << std::endl;
		std::cout << "finite difference: " << (energy1 - energy) / eps << std::endl;
		std::cout << "directional derivative: " << deriv.dot(dir) << std::endl;
		std::cout << "error: " << std::abs((energy1 - energy) / eps - deriv.dot(dir)) << std::endl;

		std::cout << std::endl << "derivative - hessian: " << std::endl;
		std::cout << "finite difference: " << (deriv1 - deriv).norm() / eps << std::endl;
		std::cout << "directional derivative: " << (hessian * dir).norm() << std::endl;
		std::cout << "error: " << ((deriv1 - deriv) / eps - hessian * dir).norm() << std::endl;
	}
	_state = backupState;
	_isPosHess = backupPosHess;
}

void WTFShell::testElasticReducedEnergy()
{
	auto backupPosHess = _isPosHess;
	_isPosHess = false;

	Eigen::VectorXd deriv;
	Eigen::SparseMatrix<double> hessian;

	double energy = elasticReducedEnergy(&deriv, &hessian);

	std::cout.precision(10);
	std::cout << std::endl << "Test total energy: " << std::endl;

	Eigen::VectorXd dir = Eigen::VectorXd::Random(deriv.size());
	dir.normalize();
	auto backupState = _state;

	int nverts = _state.basePos.rows();
	int nedges = _state.baseMesh.nEdges();

	std::cout << "verts: " << nverts << std::endl;
	std::cout << "edges: " << nedges << std::endl;
	std::cout << "energy: " << energy << std::endl;

	for (int i = 3; i < 10; i++)
	{
		double eps = std::pow(10, -i);
		for (int j = 0; j < nverts; j++)
		{
			_state.amplitude(j) = backupState.amplitude(j) + eps * dir(j);
		}
		for (int j = 0; j < nedges; j++)
		{
			_state.dphi(j) = backupState.dphi(j) + eps * dir(j + nverts);
		}

		Eigen::VectorXd deriv1;
		double energy1 = elasticReducedEnergy(&deriv1, NULL);
		//std::cout<<energy1<<std::endl;
		std::cout << std::endl << "eps: " << eps << std::endl;
		std::cout << std::endl << "energy - derivative: " << std::endl;
		std::cout << "finite difference: " << (energy1 - energy) / eps << std::endl;
		std::cout << "directional derivative: " << deriv.dot(dir) << std::endl;
		std::cout << "error: " << std::abs((energy1 - energy) / eps - deriv.dot(dir)) << std::endl;

		std::cout << std::endl << "derivative - hessian: " << std::endl;
		std::cout << "finite difference: " << (deriv1 - deriv).norm() / eps << std::endl;
		std::cout << "directional derivative: " << (hessian * dir).norm() << std::endl;
		std::cout << "error: " << ((deriv1 - deriv) / eps - hessian * dir).norm() << std::endl;
	}
	_state = backupState;
	_isPosHess = backupPosHess;
}