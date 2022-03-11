#ifndef DATACONVERSION_H 
#define DATACONVERSION_H

#include <Eigen/Dense>

void matToVec(const Eigen::MatrixXd &mat, Eigen::VectorXd &vec);
void vecToMat(const Eigen::VectorXd &vec, Eigen::MatrixXd &mat);

#endif
