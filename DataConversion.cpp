#include "DataConversion.h"

void matToVec(const Eigen::MatrixXd &mat, Eigen::VectorXd &vec)
{
     int nverts = mat.rows();
     vec.resize(3 * nverts);
     for (int i = 0; i < nverts; i++)
     {
         for (int j = 0; j < 3; j++)
             vec(3 * i + j) = mat(i,j);
     }
}

void vecToMat(const Eigen::VectorXd &vec, Eigen::MatrixXd &mat)
{
     int nverts = vec.size()/3;
     mat.resize(nverts, 3);
     for (int i = 0; i < nverts; i++)
     {
         for (int j = 0; j < 3; j++)
             mat(i, j) = vec(3 * i + j);
     }
}
