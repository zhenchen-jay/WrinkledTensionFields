#ifndef OBSTACLE_H
#define OBSTACLE_H

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <igl/per_vertex_normals.h>
#include <vector>
#include <set>

class Obstacle
{
    public : 
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    Obstacle(const Eigen::MatrixXd &obs_V, const Eigen::MatrixXi &obs_F) : V(obs_V), F(obs_F) {}
};

/*
class Obstacle{

    public :
    Obstacle(){} 
    Obstacle(const Eigen::MatrixXd obs_V, const Eigen::MatrixXi obs_F, const int numx, const int numy, const int numz) : nx(numx), ny(numy), nz(numz)
    {
        centroid.setZero();
        int nverts = obs_V.rows();

        for (int i = 0; i < nverts; i++)
            centroid += obs_V.row(i);

        centroid /= nverts;
  
        for (int i = 0; i < 3; i++)
            boxLength[i] = 1.01 * (obs_V.col(i).maxCoeff() - obs_V.col(i).minCoeff());

        //boxLength[0] = 0.55;
        gridSize[0] = double(boxLength[0]/(nx-1));
        gridSize[1] = double(boxLength[1]/(ny-1));
        gridSize[2] = double(boxLength[2]/(nz-1));
     
        
        generateGridPoints();
        generateGridSDF(obs_V, obs_F);
    } 

    Eigen::MatrixXd grid;
    Eigen::VectorXd gridvals;
    void generateGridPoints(); //generate a grid with nx x ny x nz points;
    void generateGridSDF(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
    bool isInsideGrid(const Eigen::RowVector3d &p) const;
    double gridSize[3]; 
    double boxLength[3]; 
    std::vector<double> computeSDF(const Eigen::MatrixXd &V, std::vector<std::vector<double> > &barys, std::vector<std::vector<double> > &gridSDFs) const;

    private :
        int nx; 
        int ny; 
        int nz;
        Eigen::Vector3d centroid;
};
*/

#endif
