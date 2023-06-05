
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

#include <random>
#include <iostream>
#include <vector>
#include <set>

#include <igl/cotmatrix.h>
#include <igl/boundary_loop.h>

#include "ElasticShellMaterial.h"
#include "../MeshLib/MeshConnectivity.h"
#include "../MeshLib/GeometryDerivatives.h"

double ElasticShellMaterial::cotan_v0(const Eigen::Vector3d v0, const Eigen::Vector3d v1, const Eigen::Vector3d v2)
{
    double e0 = (v2 - v1).norm(); 
    double e1 = (v2 - v0).norm(); 
    double e2 = (v0 - v1).norm(); 
    double angle0 = acos((e1 * e1 + e2 * e2 - e0 * e0) / (2 * e1 * e2));
    double cot = 1.0 / tan(angle0);
        
    return cot;
}






