#include "findCorners.h"
#include "MeshLib/MeshConnectivity.h"
#include <deque>
#include <vector>
#include <iostream>
#include <igl/boundary_loop.h>

static int labelComponents(const Eigen::MatrixXi& F, Eigen::VectorXi& labels)
{
	MeshConnectivity mesh(F);
	int nfaces = F.rows();
	labels.resize(nfaces);
	std::vector<bool> visited(nfaces);
	int label = 0;
	for (int i = 0; i < nfaces; i++)
	{
		if (visited[i])
			continue;
		std::deque<int> tovisit;
		tovisit.push_back(i);
		while (!tovisit.empty())
		{
			int next = tovisit.front();
			tovisit.pop_front();
			if (visited[next])
				continue;
			visited[next] = true;
			labels[next] = label;
			for (int j = 0; j < 3; j++)
			{
				int e = mesh.faceEdge(next, j);
				int o = mesh.faceEdgeOrientation(next, j);
				int opp = mesh.edgeFace(e, 1 - o);
				if (opp != -1 && !visited[opp])
					tovisit.push_back(opp);
			}
		}
		label++;
	}

	return label;
}

void findCorners(const Eigen::MatrixXd &V2D, const Eigen::MatrixXi& F2D, const Eigen::MatrixXd& V3D, const Eigen::MatrixXi& F3D, std::set<int>& corners)
{
	corners.clear();

	Eigen::VectorXi labels;
	int components = labelComponents(F2D, labels);
	std::cout << "2D mesh has " << components << " components." << std::endl;

	int nfaces = F3D.rows();
	int n2Dverts = V2D.rows();
	int n3Dverts = V3D.rows();
	MeshConnectivity mesh2D(F2D);
	MeshConnectivity mesh3D(F3D);

	std::vector<std::set<int> > vertlabels(n3Dverts);

	for (int i = 0; i < nfaces; i++)
	{
		int l = labels[i];
		for (int j = 0; j < 3; j++)
		{
			int v = mesh3D.faceVertex(i, j);
			vertlabels[v].insert(l);
		}
	}
	std::vector<std::vector<int> > L;
	igl::boundary_loop(F3D, L);

	for (auto& loop : L)
	{
		for (auto it : loop)
		{
			vertlabels[it].insert(-1);
		}
	}

	for (int i = 0; i < nfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int v2d = mesh2D.faceVertex(i, j);
			int v3d = mesh3D.faceVertex(i, j);
			if (vertlabels[v3d].size() > 2)
				corners.insert(v2d);
		}
	}
}

