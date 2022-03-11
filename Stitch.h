#ifndef STITCH_H
#define STITCH_H

#include <Eigen/Core>
#include "CommonFunctions.h"

void stitchMeshes(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& stitchedV, Eigen::MatrixXi& stitchedF, std::vector<Eigen::Vector3i>& bnd_edges, Eigen::VectorXi& newIndex);
void stitchMeshesWithTol(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXd &stitchedV, Eigen::MatrixXi &stitchedF, double tol);
void testStitchMeshes();
//void stitchMeshes(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXd &stitchedV, Eigen::MatrixXi &stitchedF);
#endif
  
  
 
  
 
  
  
  
  
 
  
  
 
  
  
  
  
  
  
  
  
  
 
  
  
  
  
  
 
  
  
 
  
  
  
  
  
  
