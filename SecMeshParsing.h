#pragma once
#include <map>
#include "CommonTools.h"

std::map<std::pair<int, int>, int> he2Edge(const Eigen::MatrixXi& faces);
std::map<std::pair<int, int>, int> he2Edge(const std::vector< std::vector<int>>& edgeToVert);
Eigen::VectorXd swapEdgeVec(const std::vector< std::vector<int>>& edgeToVert, const Eigen::VectorXd& edgeVec, int flag);
Eigen::VectorXd swapEdgeVec(const Eigen::MatrixXi& faces, const Eigen::VectorXd& edgeVec, int flag);
std::vector<std::vector<int>> swapEdgeIndices(const Eigen::MatrixXi& faces, const std::vector<std::vector<int>>& edgeIndices, int flag);
Eigen::MatrixXd edgeVec2FaceVec(const Mesh& mesh, Eigen::VectorXd& edgeVec);
Mesh convert2SecMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
void parseSecMesh(const Mesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F);

