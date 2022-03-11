#ifndef PENALTYENERGY_H
#define PENALTYENERGY_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include "ElasticShell/ElasticSetup.h"
#include "Obstacle.h"
#include "CommonFunctions.h"
/*
class Distance
{
 public:
  // Efficiently calculates whether or not a point p is closer to than eta distance to the plane spanned by vertices q0, q1, and q2.
  // This method does not require any floating point divisions and so is significantly faster than computing the distance itself.
  static bool vertexPlaneDistanceLessThan(const Eigen::Vector3d &p, 
					    const Eigen::Vector3d &q0, const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, double eta)
  {
	Eigen::Vector3d c = (q1-q0).cross(q2-q0);
	return c.dot(p-q0)*c.dot(p-q0) < eta*eta*c.dot(c);
  }

  // Efficiently calculates whether or not the line spanned by vertices (p0, p1) is closer than eta distance to the line spanned by vertices (q0, q1).
  // This method does not require any floating point divisions and so is significantly faster than computing the distance itself.
  static bool lineLineDistanceLessThan(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1,
					  const Eigen::Vector3d &q0, const Eigen::Vector3d &q1, double eta)
  {
	Eigen::Vector3d c = (p1-p0).cross(q1-q0);
	return c.dot(q0-p0)*c.dot(q0-p0) < eta*eta*c.dot(c);
  }

  // Computes the vector between a point p and the closest point to p on the triangle (q0, q1, q2). Also returns the barycentric coordinates of this closest point on the triangle;
  // q0bary is the barycentric coordinate of q0, etc. (The distance from p to the triangle is the norm of this vector.)
  static Eigen::Vector3d vertexFaceDistance(const Eigen::Vector3d &p, 
					    const Eigen::Vector3d &q0, const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, 
					    double &q0bary, double &q1bary, double &q2bary, Eigen::Matrix3d *dq = NULL)
  {
  Eigen::Vector3d ab = q1-q0;
  Eigen::Vector3d ac = q2-q0;
  Eigen::Vector3d ap = p-q0;

  double d1 = ab.dot(ap);
  double d2 = ac.dot(ap);


  if (dq)
      dq->setZero();

  // corner and edge cases

  if(d1 <= 0 && d2 <= 0)
    {
      q0bary = 1.0;
      q1bary = 0.0;
      q2bary = 0.0;
      return q0-p;
    }

  Eigen::Vector3d bp = p-q1;
  double d3 = ab.dot(bp);
  double d4 = ac.dot(bp);
  if(d3 >= 0 && d4 <= d3)
    {
      q0bary = 0.0;
      q1bary = 1.0;
      q2bary = 0.0;
      return q1-p;
    }

  double vc = d1*d4 - d3*d2;
  if((vc <= 0) && (d1 >= 0) && (d3 <= 0))
    {
      double v = d1 / (d1-d3);
      q0bary = 1.0 - v;
      q1bary = v;
      q2bary = 0;
      if (dq)
      {
          dq->row(0) = -ab.transpose() / (d1 - d3);
          dq->row(1) = ab.transpose() / (d1 - d3);
      }

      // notice that d1 - d3 = || q1 - q0 ||^2, and d1 is linear in q, thus hq is always 0.


      return (q0 + v*ab)-p;
    }
  
  Eigen::Vector3d cp = p-q2;
  double d5 = ab.dot(cp);
  double d6 = ac.dot(cp);
  if(d6 >= 0 && d5 <= d6)
    {
      q0bary = 0;
      q1bary = 0;
      q2bary = 1.0;
      return q2-p;
    }

  double vb = d5*d2 - d1*d6;
  if((vb <= 0) && (d2 >= 0) && (d6 <= 0))
    {
      double w = d2/(d2-d6);
      q0bary = 1-w;
      q1bary = 0;
      q2bary = w;
      if (dq)
      {
          dq->row(0) = -ac.transpose() / (d2 - d6);
          dq->row(2) = ac.transpose() / (d2 - d6);
      }

      // notice that d2 - d6 = || q2 - q0 ||^2, and d2 is linear in q, thus hq is always 0.

      return (q0 + w*ac)-p;
    }

  double va = d3*d6 - d5*d4;
  if((va <= 0) && (d4-d3 >= 0) && (d5-d6 >= 0))
    {
      double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
      q0bary = 0;
      q1bary = 1.0 - w;
      q2bary = w;
      if (dq)
      {
          dq->row(1) = (ab - ac).transpose() / (d4 - d3 + d5 - d6);
          dq->row(2) = (ac - ab).transpose() / (d4 - d3 + d5 - d6);
      }

      // notice that d4 - d6 + d5 - d3 = || q2 - q0 ||^2 + (q1 - q0)^T(q1 - q2), and d4 - d3 is linear in q, thus hq is always 0.

      return (q1 + w*(q2-q1))-p;
    }

  // face case
  double denom = 1.0 / (va + vb + vc);
  double v = vb * denom;
  double w = vc * denom;
  double u = 1.0 - v - w;
  q0bary = u;
  q1bary = v;
  q2bary = w;
  Eigen::Vector3d dv = denom * ((d5 - d1) * ac + (d2 - d6) * ab); 
  Eigen::Vector3d dw = denom * ((d1 - d3) * ac + (d4 - d2) * ab); 
  if (dq)
  {
      dq->row(0) = -(dv+dw).transpose();
      dq->row(1) = dv.transpose();   
      dq->row(2) = dw.transpose(); 
      std::cout<<dv.transpose()<<std::endl; 

       // notice that va, vb, vc is linear w.r.t. p and their sum is const w.r.t. p (just do simple calculation), thus hq is always 0. Or projection is actually a linear map. 
  }
  return (u*q0 + v*q1 + w*q2)-p;
  }

  // Computes the shotest vector between a segment (p0, p1) and segment (q0, q1). Also returns the barycentric coordinates of the closest points on both segments; p0bary is the barycentric
  // coordinate of p0, etc. (The distance between the segments is the norm of this vector).
  static Eigen::Vector3d edgeEdgeDistance(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1,
					  const Eigen::Vector3d &q0, const Eigen::Vector3d &q1,
					  double &p0bary, double &p1bary,
					  double &q0bary, double &q1bary)
  {  
  Eigen::Vector3d d1 = p1-p0;
  Eigen::Vector3d d2 = q1-q0;
  Eigen::Vector3d r = p0-q0;
  double a = d1.squaredNorm();
  double e = d2.squaredNorm();
  double f = d2.dot(r);

  double s,t;

  double c = d1.dot(r);
  double b = d1.dot(d2);
  double denom = a*e-b*b;
  if(denom != 0.0) 
    {
      s = clamp( (b*f-c*e)/denom );
    }
  else 
    {
      //parallel edges and/or degenerate edges; values of s doesn't matter
      s = 0;
    }
  double tnom = b*s + f;
  if(tnom < 0 || e == 0)
    {
      t = 0;
      if(a == 0)
	s = 0;
      else
	s = clamp(-c/a);  
    }
  else if(tnom > e)
    {
      t = 1.0;
      if(a == 0)
	s = 0;
      else
	s = clamp( (b-c)/a );
    }
  else
    t = tnom/e;	    

  Eigen::Vector3d c1 = p0 + s*d1;
  Eigen::Vector3d c2 = q0 + t*d2;

  p0bary = 1.0-s;
  p1bary = s;
  q0bary = 1.0-t;
  q1bary = t;

  return c2-c1;
  }

  // Computes the shortest distance between a triangle mesh and itself, i.e. the shortest distance between a vertex and a face that does not contain that vertex, or of an edge and another edge that
  // does not share a vertex. The mesh has verts1.size()/3 vertices, stored as consecutive triplets in the vector verts, and faces stored as vertex indices in the columns of faces. This method assumes
  // that each vertex is part of at least one triangle.
  // Distances between primitives *all* of whose vertices are in fixedVerts are ignored.
  static double meshSelfDistance(const Eigen::VectorXd &verts, const Eigen::Matrix3Xi &faces, const std::set<int> &fixedVerts);

 private:
  static double clamp(double u)
  {
    return std::min(1.0, std::max(u, 0.0));
  }
};

*/

double penaltyEnergy_Sphere(
    const Eigen::MatrixX3d &V, // Current (deformed) positions of the vertices
    double penaltyK, // penalty stiffness
    Eigen::Vector3d ctr, // radius of the ball centered at (0,0,-r)
    double r, // radius of the ball centered at (0,0,-r)
    Eigen::VectorXd *dEnergy, // If not null, the energy gradient will be written to this vector
    std::vector<Eigen::Triplet<double> > *hEnergy, // If not null, the energy Hessian will be written to this vector (in sparse matrix form)
    bool isProjHess = false
);

double springEnergy_Sphere(
    const Eigen::MatrixX3d &V, // Current (deformed) positions of the vertices
    const std::set<int> &collidingPoints, //points that collides with the sphere at least once
    double penaltyK, // penalty stiffness
    Eigen::Vector3d ctr, // radius of the ball centered at (0,0,-r)
    double r, // radius of the ball centered at (0,0,-r)
    Eigen::VectorXd *dEnergy, // If not null, the energy gradient will be written to this vector
    std::vector<Eigen::Triplet<double> > *hEnergy, // If not null, the energy Hessian will be written to this vector (in sparse matrix form)}
    bool isProjHess = false
);

double penaltyEnergy_Obstacle(
    const Eigen::MatrixXd &V, // Current (deformed) positions of the vertices
    const Eigen::MatrixXi &F, // Current (deformed) positions of the vertices
    const Eigen::MatrixXd &obs_V,
    const Eigen::MatrixXi &obs_F,
    double penaltyK, // penalty stiffness
    Eigen::VectorXd *dE, // If not null, the energy gradient will be written to this vector
    std::vector<Eigen::Triplet<double> > *hE,
    bool isProjHess = false);
//    Eigen::MatrixXd &LocalNearestP,
//    Eigen::VectorXd &LocalMinDis);

double penaltyEnergy_prescribed(
    const Eigen::MatrixXd &V, // Current (deformed) positions of the vertices
    const Eigen::MatrixXi &F,
    double penaltyK, // penalty stiffness
    Eigen::MatrixXd LocalNearestP,
    Eigen::VectorXd LocalMinDis);
//This is only C0 continuous 
//double penaltyForce_smoothedSDF(
//    const Eigen::MatrixXd &V, // Current (deformed) positions of the vertices
//    const Obstacle &obs,
//    double penaltyK, // penalty stiffness
//    Eigen::VectorXd *dE, // If not null, the energy gradient will be written to this vector
//    std::vector<Eigen::Triplet<double> > *hE,
//    bool isProjHess = false);

double penaltyForce_VertexFace(
    const Eigen::MatrixXd &cloth_start,
    const ElasticSetup &setup,
    Eigen::VectorXd *dE, // If not null, the energy gradient will be written to this vector
    std::vector<Eigen::Triplet<double> > *hE, // If not null, the energy Hessian will be written to this vector (in sparse matrix form),
    bool isProjHess);

double vertexFaceEnergy(
    const Eigen::Vector3d& p,
    const Eigen::Vector3d& q0,
    const Eigen::Vector3d& q1,
    const Eigen::Vector3d& q2,
    const double &threshold,
    const double &coeff,
    Eigen::Vector3d *deriv,
    Eigen::Matrix3d *hess
);

void testPenaltyGradient(
    const Eigen::MatrixXd &V,
    const ElasticSetup &setup);

void testVertexFaceEnergy();

#endif
