#include "PenaltyEnergy.h"
#include "../Collision/Collision.h"
#include "../Collision/Distance.h"
#include <igl/signed_distance.h>
#include <cmath>
#include <set>
#include <memory>


double penaltyEnergy_Sphere(
	const Eigen::MatrixX3d &V, // Current (deformed) positions of the vertices
	double penaltyK, // penalty stiffness
	Eigen::Vector3d ctr, // radius of the ball centered at (0,0,-r)
	double r, // radius of the ball centered at (0,0,-r)
	Eigen::VectorXd *dEnergy, // If not null, the energy gradient will be written to this vector
	std::vector<Eigen::Triplet<double> > *hEnergy, // If not null, the energy Hessian will be written to this vector (in sparse matrix form)
	bool isProjHess
)
{
	int nverts = (int)V.rows();
	double energy = 0;
	if (dEnergy)
	{
		dEnergy->resize(3 * nverts);
		dEnergy->setZero();
	}
	if (hEnergy)
		hEnergy->clear();

	for (int i = 0; i < nverts; i++)
	{
		Eigen::Vector3d diff = V.row(i).transpose() - ctr;
		if (diff.norm() < r)
		{
			energy += 0.5 * penaltyK * (r - diff.norm()) * (r - diff.norm());
			if (dEnergy)
				dEnergy->segment<3>(3 * i) -= penaltyK * (r - diff.norm()) * diff.normalized();
			if (hEnergy)
			{
				Eigen::Matrix3d localH;
				localH = -penaltyK * ((r / diff.norm() - 1.0) * Eigen::MatrixXd::Identity(3,3) - r * diff * diff.transpose() / pow(diff.norm(),3));
				if(isProjHess) 
					localH = lowRankApprox(localH);
				for (int j = 0; j < 3; j++)
				{
					for (int k = 0; k <3; k++)
					{
						hEnergy->push_back(Eigen::Triplet<double>(3 * i + j, 3 * i + k, localH(j,k)));
					}
				}
			}
		}
	}
	return energy;
}


double springEnergy_Sphere(
	const Eigen::MatrixX3d &V, // Current (deformed) positions of the vertices
	const std::set<int> &collidingPoints, //points that collides with the sphere at least once
	double penaltyK, // penalty stiffness
	Eigen::Vector3d ctr, // radius of the ball centered at (0,0,-r)
	double r,
	Eigen::VectorXd *dEnergy, // If not null, the energy gradient will be written to this vector
	std::vector<Eigen::Triplet<double> > *hEnergy, // If not null, the energy Hessian will be written to this vector (in sparse matrix form)
	bool isProjHess
)
{
	int nverts = (int)V.rows();
	double energy = 0;
	if (dEnergy)
	{
		dEnergy->resize(3 * nverts);
		dEnergy->setZero();
	}
	if (hEnergy)
		hEnergy->clear();

	for (std::set<int>::iterator it = collidingPoints.begin(); it != collidingPoints.end(); ++it)
	{
		Eigen::Vector3d diff = V.row(*it).transpose() - ctr;
			energy += 0.5 * penaltyK * (r - diff.norm()) * (r - diff.norm());
			if (dEnergy)
				dEnergy->segment<3>(3 * (*it)) -= penaltyK * (r - diff.norm()) * diff.normalized();
			if (hEnergy)
			{
				Eigen::Matrix3d localH;
				localH = -penaltyK * ((r / diff.norm() - 1.0) * Eigen::MatrixXd::Identity(3,3) - r * diff * diff.transpose() / pow(diff.norm(),3));
				if(isProjHess) 
					localH = lowRankApprox(localH);
				for (int j = 0; j < 3; j++)
				{
					for (int k = 0; k <3; k++)
					{
						hEnergy->push_back(Eigen::Triplet<double>(3 * (*it) + j, 3 * (*it) + k, localH(j,k)));
					}
				}
			}
	}
	return energy;
}

double penaltyEnergy_Obstacle(
	const Eigen::MatrixXd &V, // Current (deformed) positions of the vertices
	const Eigen::MatrixXi &F, // Current (deformed) positions of the vertices
	const Eigen::MatrixXd &obs_V,
	const Eigen::MatrixXi &obs_F,
	double penaltyK, // penalty stiffness
	Eigen::VectorXd *dE, // If not null, the energy gradient will be written to this vector
	std::vector<Eigen::Triplet<double> > *hE, // If not null, the energy Hessian will be written to this vector (in sparse matrix form),
//    Eigen::MatrixXd &LocalNearestP,
//    Eigen::VectorXd &LocalMinDis
	bool isProjHess
)
{
	 double energy = 0;
	 int nVerts = V.rows();
	 
	 if (dE)
	 {
		 dE->resize(3 * nVerts);
		 dE->setZero();
	 }
	 if (hE)
		 hE->clear();

	 if (obs_V.rows() > 0)
	 { 
		 Eigen::MatrixXd LocalNearestP, normals;
		 Eigen::VectorXi nearestF;
		 Eigen::VectorXd LocalMinDis;
		 igl::SignedDistanceType type = igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_PSEUDONORMAL;
		 
		 igl::signed_distance(V, obs_V, obs_F, type, LocalMinDis, nearestF, LocalNearestP, normals);
		 for (int j = 0; j < nVerts; j++){
			  if (LocalMinDis(j) < 0){
				   Eigen::Vector3d diff = (V.row(j)- LocalNearestP.row(j)).transpose();
				   energy += 0.5 * penaltyK * diff.squaredNorm();
				   if (dE)
					   dE->segment<3> (3 * j) += penaltyK * diff;
				   if (hE) {
					   for (int m = 0; m < 3; m++)
							hE->push_back(Eigen::Triplet<double>(3 * j + m, 3 * j + m, penaltyK));
				   }
			  }
		 }

		 // upsample faces for additional collision points
		 int nfaces = F.rows();
		 std::vector<std::vector<int> > vid(nfaces, std::vector<int>(3));  
		 //int nquads = quadPoints.size();
		 int nquads = 1;
		 std::vector<std::vector<double> > weight(nquads);
		 for (int j = 0; j < nquads; j++)
		 {
			 double u = 1/3.0;
			 double v = 1/3.0;
			 weight[j].push_back(1 - u -v);
			 weight[j].push_back(u);
			 weight[j].push_back(v);
		 }
   
		 Eigen::MatrixXd additionV(nfaces * nquads,3);
 
		 for (int i = 0; i < nfaces; i++)
		 {
			 for (int j = 0; j <3; j++)
				 vid[i][j] = F(i,j);
			 for (int j = 0; j < nquads; j++)
			 {
				 additionV.row(nquads * i + j) = weight[j][0] * V.row(vid[i][0]) + weight[j][1] * V.row(vid[i][1]) + weight[j][2] * V.row(vid[i][2]);
			 }
		}

		igl::signed_distance(additionV, obs_V, obs_F, type, LocalMinDis, nearestF, LocalNearestP, normals);
		for (int j = 0; j < additionV.rows(); j++){
			 if (LocalMinDis(j) < 0){
				 Eigen::Vector3d diff = (additionV.row(j)- LocalNearestP.row(j)).transpose();
				 energy += 0.5 * penaltyK * diff.squaredNorm();
				 if (dE)
				 {
					 for (int k = 0; k < 3; k++)
						 dE->segment<3> (3 * vid[j/nquads][k]) += penaltyK * diff * weight[j%nquads][k];
				 }
				 if (hE) {
					 Eigen::Matrix3d localH = penaltyK * Eigen::MatrixXd::Identity(3,3);
					 if(isProjHess) 
						localH = lowRankApprox(localH);
					 for (int m = 0; m < 3 ; m ++){
						 for (int l = 0; l < 3 ; l ++){
							 for (int k = 0 ; k < 3; k++)
									 hE->push_back(Eigen::Triplet<double> (3 * vid[j/nquads][m] + k , 3 * vid[j/nquads][l] + k, weight[j%nquads][m] * weight[j%nquads][l] * penaltyK));
						 }
					 } 
				 }
			 }
		 }

	}
	return energy;
}

double penaltyEnergy_prescribed(
	const Eigen::MatrixXd &V, // Current (deformed) positions of the vertices
	const Eigen::MatrixXi &F,
	double penaltyK, // penalty stiffness
	Eigen::MatrixXd LocalNearestP,
	Eigen::VectorXd LocalMinDis)
{

	 double energy = 0;
	 int nVerts = V.rows();
	 
	 for (int j = 0; j < nVerts; j++){
		 if (LocalMinDis(j) < 0){
			  Eigen::Vector3d diff = (V.row(j)- LocalNearestP.row(j)).transpose();
			  energy += 0.5 * penaltyK * diff.squaredNorm();
		 }
	 }

	  // upsample faces for additional collision points
	  int nfaces = F.rows();
	  std::vector<std::vector<int> > vid(nfaces, std::vector<int>(3));  
	  //int nquads = quadPoints.size();
	  int nquads = 1;
	  std::vector<std::vector<double> > weight(nquads);
	  for (int j = 0; j < nquads; j++)
	  {
		  double u = 1.0/3;
		  double v = 1.0/3;
		  weight[j].push_back(1 - u -v);
		  weight[j].push_back(u);
		  weight[j].push_back(v);
	  }
	  Eigen::MatrixXd additionV(nfaces * nquads,3);
 
	  for (int i = 0; i < nfaces; i++)
	  {
		  for (int j = 0; j <3; j++)
			  vid[i][j] = F(i,j);
		  for (int j = 0; j < nquads; j++)
			  additionV.row(nquads * i + j) = weight[j][0] * V.row(vid[i][0]) + weight[j][1] * V.row(vid[i][1]) + weight[j][2] * V.row(vid[i][2]);
	 }

	 for (int j = 0; j < additionV.rows(); j++){
		 if (LocalMinDis(j) < 0){
			  Eigen::Vector3d diff = (additionV.row(j)- LocalNearestP.row(j)).transpose();
			  energy += 0.5 * penaltyK * diff.squaredNorm();
		 }
	 }

	 return energy;
}

double vertexFaceEnergy(
	const Eigen::Vector3d &p,
	const Eigen::Vector3d &q0,
	const Eigen::Vector3d &q1,
	const Eigen::Vector3d &q2,
	const double &threshold,
	const double &coeff,
	Eigen::Vector3d *deriv,
	Eigen::Matrix3d *hess
)
{
	double q0bary, q1bary, q2bary;
	Eigen::Matrix3d dbarys;
	Distance::vertexFaceDistance(p, q0, q1, q2, q0bary, q1bary, q2bary, &dbarys); // notice that hbary is always 0

	Eigen::Vector3d q = q0bary * q0 + q1bary * q1 + q2bary * q2;            // projected point

	Eigen::Matrix3d dq = q0 * dbarys.row(0) + q1 * dbarys.row(1) + q2 * dbarys.row(2);
	Eigen::Matrix3d I;
	I.setIdentity();

	double energy = 0;
	if(deriv)
		deriv->setZero();
	if(hess)
		hess->setZero();
	
	if ((p - q).norm() < threshold)
	{
		Eigen::Vector3d normDeriv = ((p - q).transpose() / (p - q).norm() * (I - dq)).transpose();
		Eigen::Matrix3d normHess = (I - dq).transpose() * (p - q) * normDeriv.transpose() * ( -1.0 / (p - q).squaredNorm() ) + 1.0 / (p-q).norm() * (I - dq).transpose() * (I - dq);

		energy = 0.5 * coeff * (threshold - (p - q).norm()) * (threshold - (p - q).norm());

		if (deriv)
			*deriv = -coeff * (threshold - (p - q).norm()) * normDeriv;

		if (hess)
			*hess = coeff * normDeriv * normDeriv.transpose() - coeff * (threshold - (p - q).norm()) * normHess;
			

	   /* double d = (p - q).norm();
		energy = -coeff * (threshold - d) * (threshold - d) * std::log(d / threshold);

		if (deriv)
			*deriv = coeff * (threshold - d) * (2.0 * std::log(d / threshold) - threshold / d + 1.0) * normDeriv;

		if (hess)
		{
			*hess = coeff * (threshold - d) * (2.0 * std::log(d / threshold) - threshold / d + 1.0) * normHess + coeff * (-2.0 * std::log(d / threshold) + (threshold - d) * (threshold + 3 * d) / (d * d)) * normDeriv * normDeriv.transpose();
		}*/
			
			
	}
	return energy;
}

double penaltyForce_VertexFace(
	const Eigen::MatrixXd &V,
	const ElasticSetup &setup,
	Eigen::VectorXd *dE, // If not null, the energy gradient will be written to this vector
	std::vector<Eigen::Triplet<double> > *hE, // If not null, the energy Hessian will be written to this vector (in sparse matrix form),
	bool isProjHess
	)
{

	int nverts = V.rows();
	Eigen::VectorXd cloth_start;
	matToVec(V, cloth_start);
	if (dE)
	{
		dE->resize(3 * nverts);
		dE->setZero();
	}
	if (hE)
		hE->clear(); 

	std::vector<Eigen::MatrixXd> obs_start;
	std::vector<std::unique_ptr<AABB> > AABBs;
	for (int i = 0; i < setup.obs.size(); i++)
	{
		obs_start.push_back(setup.obs[i].V);
		AABBs.emplace_back(buildAABB(obs_start[i], setup.obs[i].V, setup.obs[i].F, setup.innerEta));
	}

	double energy = 0;
	for (int i = 0; i < nverts; i++) 
	{
		std::map<int, double>::const_iterator it = setup.clampedDOFs.find(i * 3); // only valid when all three coordinates are fixed
		if (it == setup.clampedDOFs.end())
		{
			BoundingBox sweptLine;
			for (int j = 0; j < 3; j++)
			{
				sweptLine.mins[j] = cloth_start[3 * i + j];
				sweptLine.maxs[j] = cloth_start[3 * i + j];
			}
			for (int j = 0; j < setup.obs.size(); j++) {
				if (!AABBs[j]) // obstacle with zero faces, for some reason
					continue;

				std::vector<int> hits;
				AABBs[j]->intersect(sweptLine, hits);
				//std::cout << "vertex : " << i << " number of hits : " << hits.size() << std::endl;
				for (int k = 0; k < hits.size(); k++) 
				{
					int face = hits[k];
					Eigen::RowVector3i obs_face = setup.obs[j].F.row(face);
					Eigen::Vector3d q0 = obs_start[j].row(obs_face(0)).transpose();
					Eigen::Vector3d q1 = obs_start[j].row(obs_face(1)).transpose();
					Eigen::Vector3d q2 = obs_start[j].row(obs_face(2)).transpose();
					Eigen::Vector3d p0 = cloth_start.segment<3>(3 * i);

					Eigen::Vector3d localGrad;
					Eigen::Matrix3d localH;

					double localE = vertexFaceEnergy(p0, q0, q1, q2, setup.innerEta, setup.penaltyK, dE ? &localGrad : NULL, hE ? &localH : NULL);

					energy += localE;

					if (dE)
						dE->segment<3>(3 * i) += localGrad;
					if (hE)
					{
						if(isProjHess) 
							localH = lowRankApprox(localH);
						for (int m = 0; m < 3; m++)
						{
							for (int n = 0; n <3; n++)
							{
								hE->push_back(Eigen::Triplet<double>(3 * i + m, 3 * i + n, localH(m,n)));
							}
						}
					}
				}
			}
		}
	}
	
	return energy;
}

void testPenaltyGradient(
	const Eigen::MatrixXd &V,
	const ElasticSetup &setup
)
{
	int nverts = V.rows();
	Eigen::VectorXd grad;
	Eigen::SparseMatrix<double> hess(3 * nverts, 3 * nverts);
	std::vector<Eigen::Triplet<double>> hessList;
	double E = penaltyForce_VertexFace(V, setup, &grad, &hessList, false);
	hess.setFromTriplets(hessList.begin(), hessList.end());
	
	Eigen::VectorXd dir = Eigen::VectorXd::Random(3 * nverts);
	dir.normalize();
	
	for(int i = 3; i < 10; i++)
	{
		Eigen::VectorXd grad1;
		Eigen::VectorXd grad2;
		double eps = std::pow(10, -i);
		Eigen::MatrixXd p1 = V;
		Eigen::MatrixXd p2 = V;
		
		for (int j = 0; j < nverts; j++)
		{
			p1.row(j) = V.row(j) + eps * dir.segment(3* j, 3).transpose();
			p2.row(j) = V.row(j) - eps * dir.segment(3* j, 3).transpose();
		}

		double E1 = penaltyForce_VertexFace(p1, setup, &grad1, NULL, false);
		double E2 = penaltyForce_VertexFace(p2, setup, &grad2, NULL, false);
		std::cout << "grad size : " << grad.size() << " " << grad1.size() << " " << grad2.size()<< std::endl;

		std::cout<<"eps: "<<eps<<std::endl;
		std::cout<<"gradient check error: "<<std::abs((E1 - E) / eps - grad.dot(dir))<<std::endl;
		std::cout<<"hessian check error: "<<((grad1 - grad) / eps - hess * dir).norm()<<std::endl;
		std::cout<<"hessian check error: "<<((grad2 - grad) / eps + hess * dir).norm()<<std::endl;
	}
		

}


void testVertexFaceEnergy()
{
	Eigen::Vector3d p(0.3, 0.3, 0.1);
	Eigen::Vector3d q0(0, 0, 0);
	Eigen::Vector3d q1(1, 0, 0);
	Eigen::Vector3d q2(0, 1, 0);

	double threshold  = 10;
	double coeff = 20;

	Eigen::Vector3d dir = Eigen::Vector3d::Random();
	//dir(2) = 0;
	dir.normalize();
	dir << -0.1, 0, 0;

	Eigen::Vector3d grad;
	Eigen::Matrix3d hess;

	double E = vertexFaceEnergy(p, q0, q1, q2, threshold, coeff, &grad, &hess);

	std::cout << "E: "<<E<<std::endl;
	std::cout<<"grad: "<<grad.transpose()<<std::endl;
	std::cout<<"hess: "<<std::endl<<hess<<std::endl;
	std::cout<<"dir: "<<dir.transpose()<<std::endl;

	for(int i = 3; i < 10; i++)
	{
		Eigen::Vector3d grad1;
		Eigen::Vector3d grad2;
		double eps = std::pow(10, -i);
		Eigen::Vector3d p1 = p + eps * dir;
		Eigen::Vector3d p2 = p - eps * dir;

		double E1 = vertexFaceEnergy(p1, q0, q1, q2, threshold, coeff, &grad1, NULL);
		double E2 = vertexFaceEnergy(p2, q0, q1, q2, threshold, coeff, &grad2, NULL);

		std::cout<<"eps: "<<eps<<std::endl;
		std::cout<<"gradient check error: "<<std::abs((E1 - E) / eps - grad.dot(dir))<<std::endl;
		std::cout<<"hessian check error: "<<((grad1 - grad) / eps - hess * dir).norm()<<std::endl;
		std::cout<<"hessian check error: "<<((grad2 - grad) / eps + hess * dir).norm()<<std::endl;
		//std::cout<<"hessian check error: "<<(E1 + E2 - 2 * E)/ (eps * eps) - dir.transpose() * hess * dir<<std::endl;

	}

}
