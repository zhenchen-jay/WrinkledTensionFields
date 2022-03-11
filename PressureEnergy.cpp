
#include <iostream>
#include "GeometryDerivatives.h"
#include "PressureEnergy.h"
#include "CommonFunctions.h"

#include "Timer.h"

double pressureEnergyPerface(const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& curV,
	double pressure,
	int face,
	Eigen::VectorXd* derivative,
	Eigen::MatrixXd* hessian,
	Eigen::VectorXd center, 
	bool isProjHess)
{
	if (derivative)
	{
		derivative->resize(9);
		derivative->setZero();
	}
	if (hessian)
	{
		hessian->resize(9, 9);
		hessian->setZero();
	}

	double coef = -pressure / 6.0;
	double result = 0;

    Eigen::Vector3d v0 = curV.row(F(face, 0)).transpose() - center;
    Eigen::Vector3d v1 = curV.row(F(face, 1)).transpose() - center;
    Eigen::Vector3d v2 = curV.row(F(face, 2)).transpose() - center;

	double volume = v0.cross(v1).dot(v2);
	result = coef * (v0.cross(v1).dot(v2));

	Eigen::Vector3d n[3];
	n[0] = v1.cross(v2); 
	n[1] = v2.cross(v0); 
	n[2] = v0.cross(v1); 

	if (derivative)
	{
		for (int j = 0; j < 3; j++)
		{
			derivative->segment<3>(3 * j) = coef * n[j];
		}
	}
	

	std::vector<Eigen::Matrix3d> H(3);

	H[0] = coef * crossMatrix(v0);
	H[1] = coef * crossMatrix(v1);
	H[2] = coef * crossMatrix(v2);

	if (hessian)
	{
		hessian->block<3, 3>(0, 3) = -H[2];
		hessian->block<3, 3>(0, 6) = H[1];
		hessian->block<3, 3>(3, 0) = H[2];
		hessian->block<3, 3>(3, 6) = -H[0];
		hessian->block<3, 3>(6, 0) = -H[1];
		hessian->block<3, 3>(6, 3) = H[0];

		if (isProjHess)
			*hessian = lowRankApprox(*hessian);
	}

	return result;
}

double pressureEnergy(
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& curV,
	double pressure,
	Eigen::VectorXd* dEnergy, // positions, then thetas
	std::vector<Eigen::Triplet<double> >* hessian,
	Eigen::VectorXd center,
	bool isProjHess,
	bool isParallel)
{

	int nfaces = F.rows();
	int nverts = curV.rows();

	if (dEnergy)
	{
		dEnergy->resize(3 * nverts);
		dEnergy->setZero();
	}
	if (hessian)
	{
		hessian->clear();
	}

	double result = 0;

    /*
	for (int i = 0; i < nfaces; i++)
	{
        Eigen::Vector3d v0 = curV.row(F(i, 0)).transpose() - center;
		Eigen::Vector3d v1 = curV.row(F(i, 1)).transpose() - center;
		Eigen::Vector3d v2 = curV.row(F(i, 2)).transpose() - center;
		double volume = v0.cross(v1).dot(v2);
		Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);

		Eigen::VectorXd localGrad;
		Eigen::MatrixXd localHess;

		result += pressureEnergyPerface(F, curV, pressure, i, dEnergy || hessian ? &localGrad : NULL, dEnergy || hessian ? &localHess : NULL);
		
		for (int j = 0; j < 3; j++)
		{
			if(dEnergy)
				dEnergy->segment<3>(3 * F(i, j)) += localGrad.segment<3>(3 * j);
			if (hessian)
			{
				if (isProjHess)
					localHess = lowRankApprox(localHess);

				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						hessian->push_back(Eigen::Triplet<double>(3 * F(i, j) + k, 3 * F(i, 0) + l, localHess(3 * j + k, 3 * 0 + l)));
						hessian->push_back(Eigen::Triplet<double>(3 * F(i, j) + k, 3 * F(i, 1) + l, localHess(3 * j + k, 3 * 1 + l)));
						hessian->push_back(Eigen::Triplet<double>(3 * F(i, j) + k, 3 * F(i, 2) + l, localHess(3 * j + k, 3 * 2 + l)));
					}
				}
			}
		}
	}
    */

	auto energies = std::vector<double>(nfaces);
	auto derivs = std::vector<Eigen::VectorXd>(nfaces);
	auto hesses = std::vector<Eigen::MatrixXd>(nfaces);

	if (isParallel)
	{
		auto computePressure = [&](const tbb::blocked_range<uint32_t>& range)
		{
			for (uint32_t i = range.begin(); i < range.end(); ++i)
			{
				energies[i] = pressureEnergyPerface(F, curV, pressure, i, dEnergy || hessian ? &derivs[i] : NULL, dEnergy || hessian ? &hesses[i] : NULL, center, isProjHess);
			}
		};

		tbb::blocked_range<uint32_t> rangex(0u, (uint32_t)nfaces, GRAIN_SIZE);
		tbb::parallel_for(rangex, computePressure);
	}
	else
	{
		for (int i = 0; i < nfaces; i++)
		{
			energies[i] = pressureEnergyPerface(F, curV, pressure, i, dEnergy || hessian ? &derivs[i] : NULL, dEnergy || hessian ? &hesses[i] : NULL, center, isProjHess);
		}
	}

	for (int i = 0; i < nfaces; i++)
	{
		result += energies[i];
		for (int j = 0; j < 3; j++)
		{
			if (dEnergy)
				dEnergy->segment<3>(3 * F(i, j)) += derivs[i].segment<3>(3 * j);
			if (hessian)
			{
				Eigen::MatrixXd localHess = hesses[i];

				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						hessian->push_back(Eigen::Triplet<double>(3 * F(i, j) + k, 3 * F(i, 0) + l, localHess(3 * j + k, 3 * 0 + l)));
						hessian->push_back(Eigen::Triplet<double>(3 * F(i, j) + k, 3 * F(i, 1) + l, localHess(3 * j + k, 3 * 1 + l)));
						hessian->push_back(Eigen::Triplet<double>(3 * F(i, j) + k, 3 * F(i, 2) + l, localHess(3 * j + k, 3 * 2 + l)));
					}
				}
			}
		}
	}
	return result;
}


double pressureEnergyBackup(
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& curV,
	double pressure,
	Eigen::VectorXd* dEnergy, // positions, then thetas
	std::vector<Eigen::Triplet<double> >* hessian,
	Eigen::VectorXd center)
{

	int nfaces = F.rows();
	int nverts = curV.rows();

	if (dEnergy)
	{
		dEnergy->resize(3 * nverts);
		dEnergy->setZero();
	}
	if (hessian)
	{
		hessian->clear();
	}

	double result = 0;

	for (int i = 0; i < nfaces; i++)
	{
		Eigen::Vector3d v0 = curV.row(F(i, 0)) - center;
		Eigen::Vector3d v1 = curV.row(F(i, 1)) - center;
		Eigen::Vector3d v2 = curV.row(F(i, 2)) - center;
		double volume = v0.cross(v1).dot(v2);
		result -= pressure / 6.0 * (v0.cross(v1).dot(v2));
		Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
	/*	Eigen::Vector3d n[3];
		n[0] = v1.cross(v2);
		n[1] = v2.cross(v0);
		n[2] = v0.cross(v1);*/

		if (dEnergy || hessian) {
			for (int j = 0; j < 3; j++)
			{
				if(dEnergy)
					//dEnergy->segment<3>(3 * F(i, j)) -= pressure / 6.0 * n[j];
					dEnergy->segment<3>(3 * F(i, j)) -= pressure / 6.0 * n;
				Eigen::Matrix3d H0 = -pressure / 6.0 * crossMatrix(v2 - v1);
				Eigen::Matrix3d H1 = -pressure / 6.0 * crossMatrix(v0 - v2);
				Eigen::Matrix3d H2 = -pressure / 6.0 * crossMatrix(v1 - v0);


				if (hessian) {
					for (int k = 0; k < 3; k++)
					{
						for (int l = 0; l < 3; l++)
						{
							hessian->push_back(Eigen::Triplet<double>(3 * F(i, j) + k, 3 * F(i, 0) + l, H0(k, l)));
							hessian->push_back(Eigen::Triplet<double>(3 * F(i, j) + k, 3 * F(i, 1) + l, H1(k, l)));
							hessian->push_back(Eigen::Triplet<double>(3 * F(i, j) + k, 3 * F(i, 2) + l, H2(k, l)));
						}
					}
				}
			}
		}

		
	}

	return result;
}


void testPressureEnergyPerface(
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& curPos,
	double pressure,
	int face,
	Eigen::VectorXd center
)
{
	std::cout << std::endl << "Testing pressure: " << std::endl;

	Eigen::VectorXd dp;
	Eigen::MatrixXd hp;

	double f = pressureEnergyPerface(F, curPos, pressure,face, &dp, &hp, center);

	Eigen::VectorXd dir = Eigen::VectorXd::Random(9);
	dir.normalize();
	
	for (int i = 1; i < 10; i++)
	{
		double eps = std::pow(10, -i);
		Eigen::MatrixXd disturbPos = curPos;
		Eigen::VectorXd dp1;
		for (int j = 0; j < 3; j++)
		{
			int vid = F(face, j);
			disturbPos.row(vid) += eps * dir.segment(3 * j, 3);
		}
			
		double f1 = pressureEnergyPerface(F, disturbPos, pressure, face, &dp1, NULL, center);

		std::cout << "EPS: " << eps << std::endl;
		std::cout << "energy : " << f << ", after perturbation: " << f1 << std::endl;
		std::cout << "finite difference: " << (f1 - f) / eps << ", " << "directional deriv: " << dir.dot(dp) << ", error: " << std::abs((f1 - f) / eps - dir.dot(dp)) << std::endl;
		std::cout << "finite difference: " << ((dp1 - dp) / eps).transpose()<<std::endl<< "directional deriv: " << (hp * dir).transpose()<<std::endl << "error: " << ((dp1 - dp) / eps - hp * dir).norm() << std::endl;
	}
}

void testPressureEnergy(
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& curPos,
	double pressure,
	Eigen::VectorXd center
)
{
	std::cout << std::endl << "Testing pressure: " << std::endl;

	Eigen::VectorXd dp;
	std::vector<Eigen::Triplet<double> > htp;
	Eigen::SparseMatrix<double> hp;

	double f = pressureEnergy(F, curPos, pressure, &dp, &htp, center);

	int nverts = curPos.rows();
	hp.resize(3 * nverts, 3 * nverts);
	hp.setFromTriplets(htp.begin(), htp.end());

	Eigen::VectorXd dir = Eigen::VectorXd::Random(3 * nverts);
	dir.normalize();

	for (int i = 4; i < 10; i++)
	{
		double eps = std::pow(10, -i);
		Eigen::MatrixXd disturbPos = curPos;
		Eigen::VectorXd dp1;
		for (int i = 0; i < nverts; i++)
		{
			for(int j = 0; j < 3; j++)
				disturbPos(i, j) += eps * dir(3 * i + j);
		}

		double f1 = pressureEnergy(F, disturbPos, pressure, &dp1, NULL, center);

		std::cout << "EPS: " << eps << std::endl;
		std::cout << "energy : " << f << ", after perturbation: " << f1 << std::endl;
		std::cout << "finite difference: " << (f1 - f) / eps << ", " << "directional deriv: " << dir.dot(dp) << ", error: " << std::abs((f1 - f) / eps - dir.dot(dp)) << std::endl;
		std::cout << "finite difference: " << ((dp1 - dp) / eps).norm() << "directional deriv: " << (hp * dir).norm() << "error: " << ((dp1 - dp) / eps - hp * dir).norm() << std::endl;
	}

	//Eigen::VectorXd dp0, dp1;
	//std::vector<Eigen::Triplet<double> > htp0, htp1;
	//Eigen::SparseMatrix<double> hp0, hp1;
	////int nverts = curPos.rows();

	//double p0 = pressureEnergy(F, curPos, pressure, &dp0, &htp0, center);
	//
	//hp0.resize(3 * nverts, 3 * nverts);
	//hp0.setFromTriplets(htp0.begin(), htp0.end());


	//double p1 = pressureEnergyBackup(F, curPos, pressure, &dp1, &htp1, center);


	//hp1.resize(3 * nverts, 3 * nverts);
	//hp1.setFromTriplets(htp1.begin(), htp1.end());

	//std::cout << "energy difference : " << p0 << ", " << p0 - p1 << std::endl;
	//std::cout << "gradient difference: " << dp0.norm() << ", " << (dp0 - dp1).norm() << std::endl;
	//std::cout << "hessian difference: " << hp0.norm() << ", " <<(hp0 - hp1).norm() << std::endl;




}
