#include <igl/adjacency_list.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/unique.h>
#include <igl/per_vertex_normals.h>
#include <igl/writeOBJ.h>
#include <igl/boundary_loop.h>
#include <igl/cotmatrix.h>
#include <igl/principal_curvature.h>

#include "PhiEstimate.h"
#include "KnoppelStripePatternEdgeOmega.h"
#include "MeshUpsampling.h"
#include "../SecMeshParse.h"

#include "../ComplexLoop/ComplexLoop.h"
#include "../ComplexLoop/ComplexLoopZuenko.h"
#include "../MeshLib/GeometryDerivatives.h"
#include "../MeshLib/MeshGeometry.h"
#include "../findCorners.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <iterator>

static const double PI = 3.141592653589;

struct ClusterInfo
{
	std::vector<std::vector<std::pair<int, int> > > clusters;
	std::vector<int> jump;
};


// computes the integer multiple of 2Pi jump between the values of cutFunction on edgeFace(i,0) and edgeFace(i,1)
static void computeEdgeJumps(int uncutVerts, const Eigen::MatrixXi& uncutF, const std::set<int>& problemFaces, 
	const Eigen::VectorXd &soupFunction, Eigen::VectorXi &jumps,
	std::map<int, ClusterInfo> &clusterjumps)
{
	MeshConnectivity uncutMesh(uncutF);

	int nedges = uncutMesh.nEdges();
	jumps.resize(nedges);
	jumps.setZero();
	for (int i = 0; i < nedges; i++)
	{
		int face1 = uncutMesh.edgeFace(i, 0);
		int face2 = uncutMesh.edgeFace(i, 1);
		if (face1 == -1 || face2 == -1)
			continue;
		if (problemFaces.count(face1) || problemFaces.count(face2))
			continue;
		int vert[2];
		vert[0] = uncutMesh.edgeVertex(i, 0);
		vert[1] = uncutMesh.edgeVertex(i, 1);
		double vert1val[2];
		double vert2val[2];
		double jump[2];
		for (int k = 0; k < 2; k++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (uncutF(face1, j) == vert[k])
					vert1val[k] = soupFunction[3 * face1 + j];
				if (uncutF(face2, j) == vert[k])
					vert2val[k] = soupFunction[3 * face2 + j];

			}
			jump[k] = (vert2val[k] - vert1val[k]) / 2.0 / PI;			
		}
		
		int ijmp1 = std::round(jump[0]);
		int ijmp2 = std::round(jump[1]);
		int avjmp = (ijmp1 + ijmp2) / 2;
		jumps[i] = avjmp;
		if (ijmp1 != ijmp2 || std::fabs(avjmp - jumps[i]) > 1e-6)
		{
			std::cout << "Bad jump: " << i << " " << jump[0] << " " << jump[1] << " " << face1 << " " << face2 << " " << vert[0] << " " << vert1val[0] << " " << vert2val[0] << " " << vert[1] << " " << vert1val[1] << " " << vert2val[1] << " " << problemFaces.count(face1) << " " << problemFaces.count(face2) << std::endl;
		}		
	}

	clusterjumps.clear();
	std::vector<std::vector<std::pair<int, int> > > vertexFaces(uncutVerts);
	for (int i = 0; i < uncutMesh.nFaces(); i++)
	{
		if (!problemFaces.count(i))
		{
			for (int j = 0; j < 3; j++)
			{
				vertexFaces[uncutF(i, j)].push_back({ i, j});
			}
		}
	}

	for (int i = 0; i < uncutVerts; i++)
	{
		int ncfaces = vertexFaces[i].size();
		std::vector<bool> visited(ncfaces);
		std::map<int, int> face2id;
		for (int j = 0; j < ncfaces; j++)
		{
			face2id[vertexFaces[i][j].first] = j;
		}
		std::vector < std::vector<std::pair<int, int> > > clusters;
		for (int j = 0; j < ncfaces; j++)
		{
			if (visited[j])
				continue;
			std::vector<std::pair<int, int> > newcluster;
			std::deque<int> tovisit;
			tovisit.push_back(j);
			while (!tovisit.empty())
			{
				int next = tovisit.front();
				tovisit.pop_front();
				if (visited[next])
					continue;
				visited[next] = true;
				newcluster.push_back(vertexFaces[i][next]);
				for (int k = 0; k < 3; k++)
				{
					int edge = uncutMesh.faceEdge(vertexFaces[i][next].first, k);
					int orient = uncutMesh.faceEdgeOrientation(vertexFaces[i][next].first, k);
					int opp = uncutMesh.edgeFace(edge, 1 - orient);
					if (opp != -1 && !problemFaces.count(opp))
					{
						auto it = face2id.find(opp);
						if (it != face2id.end())
						{
							if (!visited[it->second])
							{
								tovisit.push_back(it->second);
							}
						}
					}
				}
			}
			clusters.push_back(newcluster);
		}
		if (clusters.size() > 1)
		{
			ClusterInfo ci;
			ci.clusters = clusters;
			ci.jump.resize(clusters.size());
			ci.jump[0] = 0;
			for (int i = 1; i < clusters.size(); i++)
			{
				double vf1val = soupFunction[3 * ci.clusters[0][0].first + ci.clusters[0][0].second];
				double vf2val = soupFunction[3 * ci.clusters[i][0].first + ci.clusters[i][0].second];
				double jump = (vf2val - vf1val) / 2.0 / PI;
				ci.jump[i] = std::round(jump);
				if (std::fabs(jump - ci.jump[i]) > 1e-6)
				{
					std::cout << "Bad cluster jump: " << i << std::endl;
				}
			}
			clusterjumps[i] = ci;
		}
	}
}

void neighborhoodSoupJumps(const MeshConnectivity& mesh,
	const std::set<int>& problemFaces,
	const std::vector<std::vector<int> >& vertEdges,
	const Eigen::VectorXi& periodJumps,
	int face, int vertidx,
	std::map<int, int>& soupjumps,
	const std::map<int, ClusterInfo> &clusterInfo
	)
{
	soupjumps.clear();

	std::map<int, std::vector<std::pair<int, int> > > graph;
	for (auto it : vertEdges[mesh.faceVertex(face, vertidx)])
	{
		int face1 = mesh.edgeFace(it, 0);
		int face2 = mesh.edgeFace(it, 1);
		if (face1 == -1 || face2 == -1)
			continue;
		if (problemFaces.count(face1) || problemFaces.count(face2))
			continue;
		graph[face1].push_back({ face2, periodJumps[it] });
		graph[face2].push_back({ face1, -periodJumps[it] });
	}
	struct Visit
	{
		int face;
		int totaljump;
	};
	std::deque<Visit> q;
	q.push_back({ face, 0 });
	std::map<int, bool> visited;
	while (!q.empty())
	{
		Visit next = q.front();
		q.pop_front();
		if (visited[next.face])
			continue;
		visited[next.face] = true;
		soupjumps[next.face] = next.totaljump;
		for (auto nb : graph[next.face])
		{
			if (!visited[nb.first])
			{
				q.push_back({ nb.first, next.totaljump + nb.second });
			}
		}
	}

	auto it = clusterInfo.find(mesh.faceVertex(face, vertidx));
	if (it != clusterInfo.end())
	{
		int curcluster = -1;
		int basejump = 0;
		int nclusters = it->second.clusters.size();
		for (int i = 0; i < nclusters; i++)
		{
			auto sjit = soupjumps.find(it->second.clusters[i][0].first);
			if (sjit != soupjumps.end())
			{
				curcluster = i;
				basejump = sjit->second - it->second.jump[i];
			}
		}
		
		for (int i = 0; i < nclusters; i++)
		{
			if (i == curcluster)
				continue;

			int seedface = it->second.clusters[i][0].first;
						
			std::deque<Visit> q;
			q.push_back({ seedface, basejump + it->second.jump[i]});
			std::map<int, bool> visited;
			while (!q.empty())
			{
				Visit next = q.front();
				q.pop_front();
				if (visited[next.face])
					continue;
				visited[next.face] = true;
				soupjumps[next.face] = next.totaljump;
				for (auto nb : graph[next.face])
				{
					if (!visited[nb.first])
					{
						q.push_back({ nb.first, next.totaljump + nb.second });
					}
				}
			}
		}
	}
}


static void splitIsolatedProblemFaces(
	const Eigen::MatrixXd &uncutV, 
	const Eigen::MatrixXi &uncutF, 
	const std::set<int> &uncutProblemFaces, 
	const Eigen::VectorXd &uncutsoupamp, 
	const Eigen::VectorXd &uncutsoupphi,
	Eigen::MatrixXd &splitV, 
	Eigen::MatrixXi &splitF,
	std::set<int> &splitProblemFaces, 
	Eigen::VectorXd &splitSoupamp, 
	Eigen::VectorXd &splitSoupphi)
{
	int nverts = uncutV.rows();
	std::vector<bool> problemVert(nverts);
	std::vector<bool> transitionVert(nverts);
	MeshConnectivity mesh(uncutF);

	for (int i = 0; i < mesh.nEdges(); i++)
	{
		int face1 = mesh.edgeFace(i, 0);
		int face2 = mesh.edgeFace(i, 1);
		
		int numproblems = 0;
		int numok = 0;
		if (face1 != -1 && uncutProblemFaces.count(face1)) numproblems++;
		else if (face1 != -1) numok++;
		if (face2 != -1 && uncutProblemFaces.count(face2)) numproblems++;
		else if (face2 != -1) numok++;

		int vert1 = mesh.edgeVertex(i, 0);
		int vert2 = mesh.edgeVertex(i, 1);

		if (numproblems > 0)
		{
			problemVert[vert1] = true;
			problemVert[vert2] = true;
			if (numok > 0)
			{
				transitionVert[vert1] = true;
				transitionVert[vert2] = true;
			}
		}
	}

	std::set<int> tosplit;

	int nfaces = uncutF.rows();
	for (int i = 0; i < nfaces; i++)
	{
		if (!uncutProblemFaces.count(i))
			continue;
		bool ok = false;
		for (int j = 0; j < 3; j++)
		{
			if (problemVert[uncutF(i, j)] && !transitionVert[uncutF(i, j)])
				ok = true;
		}
		if (!ok) tosplit.insert(i);
	}

	int newfaces = nfaces + 2 * tosplit.size();
	splitF.resize(newfaces, 3);
	splitSoupphi.resize(3 * newfaces);
	splitSoupamp.resize(3 * newfaces);
	int newverts = nverts + tosplit.size();
	splitV.resize(newverts, 3);

	for (int i = 0; i < uncutV.rows(); i++)
		splitV.row(i) = uncutV.row(i);

	int currow = 0;
	int nextvert = nverts;
	for (int i = 0; i < nfaces; i++)
	{
		if (tosplit.count(i))
		{
			Eigen::Vector3d centroid(0, 0, 0);
			for (int j = 0; j < 3; j++)
				centroid += uncutV.row(uncutF(i, j)).transpose();
			centroid /= 3.0;

			splitV.row(nextvert) = centroid.transpose();

			for (int j = 0; j < 3; j++)
			{
				int v1 = j;
				int v2 = (j + 1) % 3;
				splitF(currow, 0) = uncutF(i,v1);
				splitF(currow, 1) = uncutF(i,v2);
				splitF(currow, 2) = nextvert;

				splitSoupamp[3 * currow + 0] = uncutsoupamp[3 * i + v1];
				splitSoupamp[3 * currow + 1] = uncutsoupamp[3 * i + v2];
				splitSoupamp[3 * currow + 2] = 0;

				splitSoupphi[3 * currow + 0] = uncutsoupphi[3 * i + v1];
				splitSoupphi[3 * currow + 1] = uncutsoupphi[3 * i + v2];
				splitSoupphi[3 * currow + 2] = 0;


				splitProblemFaces.insert(currow);

				currow++;
			}

			nextvert++;
		}
		else
		{
			splitF.row(currow) = uncutF.row(i);
			for (int j = 0; j < 3; j++)
			{
				splitSoupphi[3 * currow + j] = uncutsoupphi[3 * i + j];
				splitSoupamp[3 * currow + j] = uncutsoupamp[3 * i + j];
			}
			if (uncutProblemFaces.count(i))
				splitProblemFaces.insert(currow);
			currow++;
		}
	}
}

static void loopUncut(
	int uncutVerts,
	const Eigen::MatrixXi& uncutF,	
	const std::set<int> &cornerVerts,
	Eigen::SparseMatrix<double>& uncutS,
	Eigen::MatrixXi& newF
)
{
	MeshConnectivity mesh(uncutF);
	int nfaces = mesh.nFaces();
	int nedges = mesh.nEdges();

	std::vector<std::set<int> > vertexNeighbors(uncutVerts); // all neighbors
	std::vector<bool> boundaryVertex(uncutVerts); // on the boundary?
	std::vector<std::set<int> > boundaryNeighbors(uncutVerts); // neighbors on the boundary
	for (int i = 0; i < nfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{			
			vertexNeighbors[uncutF(i, j)].insert(uncutF(i, (j + 1) % 3));
			vertexNeighbors[uncutF(i, j)].insert(uncutF(i, (j + 2) % 3));			
		}
	}

	for (int i = 0; i < nedges; i++)
	{
		int face1 = mesh.edgeFace(i, 0);
		int face2 = mesh.edgeFace(i, 1);
		int vert1 = mesh.edgeVertex(i, 0);
		int vert2 = mesh.edgeVertex(i, 1);
		if (face1 == -1 || face2 == -1)
		{
			boundaryVertex[vert1] = true;
			boundaryVertex[vert2] = true;
			boundaryNeighbors[vert1].insert(vert2);
			boundaryNeighbors[vert2].insert(vert1);
		}
	}

	// Step 2: make newF
	int newfaces = 4 * nfaces;
	int newverts = uncutVerts + nedges;
	
	newF.resize(newfaces, 3);

	for (int i = 0; i < nfaces; i++)
	{
		// the central triangle
		for (int j = 0; j < 3; j++)
		{
			newF(4 * i, j) = uncutVerts + mesh.faceEdge(i, j);
		}

		// the three corner triangles
		// vertex i, edge (i+2), edge (i+1)
		for (int j = 0; j < 3; j++)
		{
			newF(4 * i + j + 1, 0) = mesh.faceVertex(i, j);
			newF(4 * i + j + 1, 1) = uncutVerts + mesh.faceEdge(i, (j + 2) % 3);
			newF(4 * i + j + 1, 2) = uncutVerts + mesh.faceEdge(i, (j + 1) % 3);

		}
	}


	// Step 3: make the uncutS stencil
	std::vector<Eigen::Triplet<double> > uncutScoeffs;
	// the old vertices
	for (int i = 0; i < uncutVerts; i++)
	{
		if (boundaryVertex[i])
		{
			if (cornerVerts.count(i))
			{
				uncutScoeffs.push_back({ i,i,1.0 });
			}
			else
			{
				int valence = boundaryNeighbors[i].size();
				assert(valence >= 2);
				double beta = 1.0 / 4.0 / double(valence);
				uncutScoeffs.push_back({ i,i,3.0 / 4.0 });
				for (auto it : boundaryNeighbors[i])
				{
					uncutScoeffs.push_back({ i, it, beta });
				}
			}
		}
		else
		{
			int valence = vertexNeighbors[i].size();
			assert(valence >= 3);
			double beta = 0;
			if (valence == 3)
			{
				beta = 3.0 / 16.0;
			}
			else
			{
				beta = 3.0 / 8.0 / double(valence);
			}

			uncutScoeffs.push_back({ i, i, 1.0 - double(valence) * beta });
			for (auto it : vertexNeighbors[i])
			{
				uncutScoeffs.push_back({ i, it, beta });
			}
		}
	}
	// the new vertices
	for (int i = 0; i < nedges; i++)
	{
		int face1 = mesh.edgeFace(i, 0);
		int face2 = mesh.edgeFace(i, 1);
		int vert1 = mesh.edgeVertex(i, 0);
		int vert2 = mesh.edgeVertex(i, 1);
		if (face1 == -1 || face2 == -1)
		{
			uncutScoeffs.push_back({ uncutVerts + i, vert1, 0.5 });
			uncutScoeffs.push_back({ uncutVerts + i, vert2, 0.5 });
		}
		else
		{
			uncutScoeffs.push_back({ uncutVerts + i, vert1, 3.0/8.0 });
			uncutScoeffs.push_back({ uncutVerts + i, vert2, 3.0/8.0 });
			int oppvert1 = mesh.edgeOppositeVertex(i, 0);
			int oppvert2 = mesh.edgeOppositeVertex(i, 1);
			uncutScoeffs.push_back({ uncutVerts + i, oppvert1, 1.0 / 8.0 });
			uncutScoeffs.push_back({ uncutVerts + i, oppvert2, 1.0 / 8.0 });
		}
	}
	uncutS.resize(newverts, uncutVerts);
	uncutS.setFromTriplets(uncutScoeffs.begin(), uncutScoeffs.end());
}

// Applies Loop subdivision to a given mesh uncutF.
// Returns three stencils:
//  - uncutS: the subdivision stencil for uncutF (ordinary Loop subdivision)
//  - soupS: the subdivision stencil for a triangle soup with the same connectivity as uncutF (only makes sense if the function at coincident
//           vertices actually agree).
//  - periodicS, periodicb: affine stencil for a periodic function on the triangle soup whose period jumps on edges is given by periodJumps.
//                                  Extends the function along edges into the "problem faces" (and otherwise treats them as if they didn't exist)
static void loopWithSeams(
	int uncutVerts,
	const Eigen::MatrixXi& uncutF,
	const Eigen::VectorXi& periodJumps,
	const std::map<int, ClusterInfo>& clusterJumps,
	const std::set<int>& problemFaces,
	const std::set<int>& cornerVerts,	
	Eigen::SparseMatrix<double>& uncutS,
	Eigen::SparseMatrix<double>& soupS,
	Eigen::SparseMatrix<double>& periodicS,
	Eigen::VectorXd& periodicb,
	Eigen::MatrixXi& newF,
	std::set<int>& newProblemFaces)
{
	MeshConnectivity mesh(uncutF);

	// Step 1: collect neighbors

	int nfaces = uncutF.rows();
	int nedges = mesh.nEdges();

	std::vector<std::set<int> > vertexNeighbors(uncutVerts); // all neighbors
	std::vector<bool> boundaryVertex(uncutVerts); // on the boundary?
	std::vector<std::set<int> > boundaryNeighbors(uncutVerts); // neighbors on the boundary
	std::vector<bool> problemVertex(uncutVerts); // touching a problem face?
	std::vector<bool> transitionVertex(uncutVerts); // on the edge between a regular and problem face?
	std::vector<bool> protectedEdges(nedges);
	/*
	for (int i = 0; i < uncutF.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (uncutF(i, j) < 0 || uncutF(i, j) >= uncutVerts)
				exit(0);
		}
	}
	*/

	for (int i = 0; i < nfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			vertexNeighbors[uncutF(i, j)].insert(uncutF(i, (j + 1) % 3));
			vertexNeighbors[uncutF(i, j)].insert(uncutF(i, (j + 2) % 3));
		}
	}

	for (int i = 0; i < nedges; i++)
	{
		int face1 = mesh.edgeFace(i, 0);
		int face2 = mesh.edgeFace(i, 1);
		int vert1 = mesh.edgeVertex(i, 0);
		int vert2 = mesh.edgeVertex(i, 1);
		if (face1 == -1 || face2 == -1)
		{
			boundaryVertex[vert1] = true;
			boundaryVertex[vert2] = true;
			boundaryNeighbors[vert1].insert(vert2);
			boundaryNeighbors[vert2].insert(vert1);
		}
		else
		{
			int problemfaces = 0;
			if (problemFaces.count(face1)) problemfaces++;
			if (problemFaces.count(face2)) problemfaces++;
			if (problemfaces == 1)
			{
				transitionVertex[vert1] = true;
				transitionVertex[vert2] = true;
			}
			if (problemfaces > 0)
			{
				problemVertex[vert1] = true;
				problemVertex[vert2] = true;
			}
		}
	}

	for (int i = 0; i < nedges; i++)
	{
		int face1 = mesh.edgeFace(i, 0);
		int face2 = mesh.edgeFace(i, 1);

		int vert1 = mesh.edgeVertex(i, 0);
		int vert2 = mesh.edgeVertex(i, 1);
		if (face1 != -1 && face2 != -1
			&& problemFaces.count(face1) && problemFaces.count(face2)
			&& transitionVertex[vert1] && transitionVertex[vert2])
		{
			// edge is chord through the problem region
			protectedEdges[i] = true;
		}
		if (((face1 == -1 && problemFaces.count(face2)) || (face2 == -1 && problemFaces.count(face1)))
			&& transitionVertex[vert1] && transitionVertex[vert2])
		{
			protectedEdges[i] = true;
		}
	}

	loopUncut(uncutVerts, uncutF, cornerVerts, uncutS, newF);

	int newfaces = newF.rows();
	int newverts = uncutVerts + nedges;

	std::vector<int> soupNewToOld(3 * newfaces + newverts, -1);
	for (int i = 0; i < nfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			soupNewToOld[3 * (4 * i + j + 1) + 0] = 3 * i + j;
		}
	}

	// Step 4: compute soup neighbor data	
	std::vector<std::vector<std::pair<int, int> > > soupNeighbors(uncutVerts); // neighbors in the triangle soup	
	std::vector<std::vector<std::pair<int, int> > > soupOrangeNeighbors(uncutVerts); // extrapolate into the problem region
	std::vector<std::vector<std::pair<int, int> > > soupBoundaryNeighbors(uncutVerts); // neighbors on the boundary or problem face in the soup
	std::vector<std::vector<std::pair<int, int> > > soupOrangeBoundaryNeighbors(uncutVerts); // extrapolate into the problem region
	std::map<int, std::pair<int, int> > nonProblemCopy;

	for (int i = 0; i < nfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{		
			if (!problemFaces.count(i))
			{
				nonProblemCopy[uncutF(i, j)] = { i,j };
			}
			int edge = mesh.faceEdge(i, j);

			soupNeighbors[uncutF(i, (j + 1) % 3)].push_back({ i,(j + 2) % 3 });
			soupNeighbors[uncutF(i, (j + 2) % 3)].push_back({ i,(j + 1) % 3 });
			
			int face1 = mesh.edgeFace(edge, 0);
			int face2 = mesh.edgeFace(edge, 1);

			if ((face1 == -1 || !problemFaces.count(face1))
				&& (face2 == -1 || !problemFaces.count(face2)))
			{
				soupOrangeNeighbors[uncutF(i, (j + 1) % 3)].push_back({ i,(j + 2) % 3 });
				soupOrangeNeighbors[uncutF(i, (j + 2) % 3)].push_back({ i,(j + 1) % 3 });
			}
			else
			{
				soupOrangeNeighbors[uncutF(i, (j + 1) % 3)].push_back({ i,(j + 1) % 3 });
				soupOrangeNeighbors[uncutF(i, (j + 2) % 3)].push_back({ i,(j + 2) % 3 });
			}

			int orient = mesh.faceEdgeOrientation(i, j);
			int oppface = mesh.edgeFace(edge, 1 - orient);
			if (oppface == -1)
			{				
				soupBoundaryNeighbors[uncutF(i, (j + 1) % 3)].push_back({ i, (j + 2) % 3 });
				soupBoundaryNeighbors[uncutF(i, (j + 2) % 3)].push_back({ i, (j + 1) % 3 });

				if (!problemFaces.count(i))
				{
					soupOrangeBoundaryNeighbors[uncutF(i, (j + 1) % 3)].push_back({ i, (j + 2) % 3 });
					soupOrangeBoundaryNeighbors[uncutF(i, (j + 2) % 3)].push_back({ i, (j + 1) % 3 });
				}
				else
				{
					soupOrangeBoundaryNeighbors[uncutF(i, (j + 1) % 3)].push_back({ i, (j + 1) % 3 });
					soupOrangeBoundaryNeighbors[uncutF(i, (j + 2) % 3)].push_back({ i, (j + 2) % 3 });
				}
			}
		}
	}
	
	// Step 5: make the soup stencil	
	std::vector<Eigen::Triplet<double> > soupScoeffs;

	for(int faceitr = 0; faceitr < newfaces; faceitr++)
	{
		for (int vertitr = 0; vertitr < 3; vertitr++)
		{
			int vertid = newF(faceitr, vertitr);
			if (vertid < uncutVerts)
			{
				if(boundaryVertex[vertid])
				{
					if (cornerVerts.count(vertid))
					{
						soupScoeffs.push_back({ 3 * faceitr + vertitr, soupNewToOld[3 * faceitr + vertitr], 1.0 });
					}
					else
					{
						int valence = soupBoundaryNeighbors[vertid].size();
						assert(valence >= 2);
						double beta = 1.0 / 4.0 / double(valence);
						double alpha = 3.0 / 4.0;

						soupScoeffs.push_back({ 3 * faceitr + vertitr, soupNewToOld[3 * faceitr + vertitr], alpha });

						for (auto it : soupBoundaryNeighbors[vertid])
						{
							soupScoeffs.push_back({ 3 * faceitr + vertitr, 3 * it.first + it.second, beta });
						}
					}
				}				
				else
				{
					int valence = soupNeighbors[vertid].size();
					assert(valence >= 6);
					double beta = 0;
					if (valence == 6)
					{
						beta = 3.0 / 32.0;
					}
					else
					{
						beta = 3.0 / 8.0 / double(valence);
					}
					double alpha = (1.0 - double(valence) * beta);
					
					soupScoeffs.push_back({ 3 * faceitr + vertitr, soupNewToOld[3 * faceitr + vertitr], alpha });
					
					for (auto it : soupNeighbors[vertid])
					{
						soupScoeffs.push_back({ 3 * faceitr + vertitr, 3 * it.first + it.second, beta });
					}
				}
			}
			else
			{
				// the new vertices
				int edgeitr = vertid - uncutVerts;

				int face1 = mesh.edgeFace(edgeitr, 0);
				int face2 = mesh.edgeFace(edgeitr, 1);
				int vert1 = mesh.edgeVertex(edgeitr, 0);
				int vert2 = mesh.edgeVertex(edgeitr, 1);
				if (face1 == -1 || face2 == -1 )
				{
					int okface = face1;
					if (okface == -1)
						okface = face2;
					for (int i = 0; i < 3; i++)
					{
						if (uncutF(okface, i) == vert1 || uncutF(okface, i) == vert2)
						{
							soupScoeffs.push_back({ 3 * faceitr + vertitr, 3 * okface + i, 0.5 });
						}
					}
				}
				else
				{
					int oppvert1 = mesh.edgeOppositeVertex(edgeitr, 0);
					int oppvert2 = mesh.edgeOppositeVertex(edgeitr, 1);
					for (int i = 0; i < 3; i++)
					{
						if (uncutF(face1, i) == oppvert1)
						{
							soupScoeffs.push_back({ 3 * faceitr + vertitr, 3 * face1 + i, 1.0 / 8.0 });
						}
						else
						{
							soupScoeffs.push_back({ 3 * faceitr + vertitr, 3 * face1 + i, 3.0 / 16.0 });
						}

						if (uncutF(face2, i) == oppvert2)
						{
							soupScoeffs.push_back({ 3 * faceitr + vertitr, 3 * face2 + i, 1.0 / 8.0 });
						}
						else
						{
							soupScoeffs.push_back({ 3 * faceitr + vertitr, 3 * face2 + i, 3.0 / 16.0 });
						}
					}
				}
			}
		}
	}
	soupS.resize(3 * newfaces, 3 * nfaces);
	soupS.setFromTriplets(soupScoeffs.begin(), soupScoeffs.end());

	// Step 6: convert the period jumps to (face, face) form
	std::map<std::pair<int, int>, int> faceFaceJumps;
	for (int i=0; i<nedges; i++)
	{
		int face0 = mesh.edgeFace(i, 0);
		int face1 = mesh.edgeFace(i, 1);
		if (face0 == -1 || face1 == -1)
			continue;
		faceFaceJumps[{face0, face1}] = periodJumps[i];
		faceFaceJumps[{face1, face0}] = -periodJumps[i];
	}

	// Step 8: compute new problem faces
	newProblemFaces.clear();
	for (auto it : problemFaces)
	{
		bool problemchild[3];
		for (int i = 0; i < 3; i++)
			problemchild[i] = false;

		for (int i = 0; i < 3; i++)
		{
			if (!transitionVertex[uncutF(it, i)])
			{
				problemchild[i] = true;
			}
			if (protectedEdges[mesh.faceEdge(it, i)])
			{
				problemchild[(i + 1) % 3] = true;
				problemchild[(i + 2) % 3] = true;
			}
		}
		int problems = 0;
		for (int i = 0; i < 3; i++)
		{
			if (problemchild[i])
			{
				problems++;
				newProblemFaces.insert(4 * it + 1 + i);
			}
		}
		if (problems != 1)
		{
			newProblemFaces.insert(4 * it);
		}
	}

	// Step 9: compute vertex edge neighors 
	std::vector<std::vector<int> > vertexEdges(uncutVerts);
	for (int i = 0; i < nedges; i++)
	{
		vertexEdges[mesh.edgeVertex(i, 0)].push_back(i);
		vertexEdges[mesh.edgeVertex(i, 1)].push_back(i);
	}	

	// Step 11: build periodic stencil

	std::vector<Eigen::Triplet<double> > periodicScoeffs;
	periodicb.resize(3 * newfaces);
	periodicb.setZero();

	for(int faceitr = 0; faceitr < newfaces; faceitr++)
	{
		if (newProblemFaces.count(faceitr))
			continue;
		
		for (int vertitr = 0; vertitr < 3; vertitr++)
		{
			int vertid = newF(faceitr, vertitr);
			if (vertid < uncutVerts)
			{
				// the old vertices

				int oldfaceid = soupNewToOld[3 * faceitr + vertitr] / 3;
				int oldvertidx = soupNewToOld[3 * faceitr + vertitr] % 3;

				if (problemFaces.count(oldfaceid))
				{
					// old vertex on a previous problem-face
					// use an arbitrary non-problem neighbor for the calculation
					oldfaceid = nonProblemCopy[vertid].first;
					oldvertidx = nonProblemCopy[vertid].second;
					vertid = uncutF(oldfaceid, oldvertidx);
				}

				std::map<int, int> soupjumps;
				neighborhoodSoupJumps(mesh, problemFaces, vertexEdges, periodJumps, oldfaceid, oldvertidx, soupjumps, clusterJumps);

				if (boundaryVertex[vertid])
				{
					if (cornerVerts.count(vertid))
					{
						periodicScoeffs.push_back({ 3 * faceitr + vertitr, 3 * oldfaceid + oldvertidx, 1.0 });
					}
					else
					{
						int valence = soupOrangeBoundaryNeighbors[vertid].size();
						assert(valence >= 2);
						double beta = 1.0 / 4.0 / double(valence);
						double alpha = 3.0 / 4.0;

						periodicScoeffs.push_back({ 3 * faceitr + vertitr, 3 * oldfaceid + oldvertidx, alpha });

						for (auto it : soupOrangeBoundaryNeighbors[vertid])
						{
							int nbface = it.first;
							int nbvert = it.second;
							if (problemFaces.count(nbface))
							{
								nbface = nonProblemCopy[uncutF(nbface, nbvert)].first;
								nbvert = nonProblemCopy[uncutF(nbface, nbvert)].second;
							}
							periodicScoeffs.push_back({ 3 * faceitr + vertitr, 3 * nbface + nbvert, beta });
							periodicb[3 * faceitr + vertitr] += -double(soupjumps[nbface]) * 2.0 * PI * beta;
						}
					}
				}
				else
				{
					int valence = soupOrangeNeighbors[vertid].size();
					assert(valence >= 6);
					double beta = 0;
					if (valence == 6)
					{
						beta = 3.0 / 32.0;
					}
					else
					{
						beta = 3.0 / 8.0 / double(valence);
					}
					double alpha = (1.0 - double(valence) * beta);

					periodicScoeffs.push_back({ 3 * faceitr + vertitr, 3*oldfaceid + oldvertidx, alpha });

					for (auto it : soupOrangeNeighbors[vertid])
					{
						int nbface = it.first;
						int nbvert = it.second;
						if (problemFaces.count(nbface))
						{
							nbface = nonProblemCopy[uncutF(nbface, nbvert)].first;
							nbvert = nonProblemCopy[uncutF(nbface, nbvert)].second;
						}
						periodicScoeffs.push_back({ 3 * faceitr + vertitr, 3 * nbface + nbvert, beta });
						periodicb[3 * faceitr + vertitr] += -double(soupjumps[nbface]) * 2.0 * PI * beta;
					}
				}
			}			
			else
			{
				// the new vertices
				int edgeitr = vertid - uncutVerts;
				int oldfaceid = faceitr / 4;

				if (problemFaces.count(oldfaceid))
				{
					// new vertex on a previously problem-face

					int corner = faceitr % 4 - 1;
					if (corner == -1)
					{
						// central triangle

						int oppface = -1;
						int oppedge = -1;
						int oppvert = -1;
						for (int i = 0; i < 3; i++)
						{
							int edge = mesh.faceEdge(oldfaceid, i);
							int edgeorient = mesh.faceEdgeOrientation(oldfaceid, i);
							int candopp = mesh.edgeFace(edge, 1 - edgeorient);
							if (candopp != -1 && !problemFaces.count(candopp))
							{
								oppface = candopp;
								oppedge = edge;
								oppvert = mesh.edgeOppositeVertex(edge, 1 - edgeorient);
							}
						}
						assert(oppface != -1);
						if (edgeitr == oppedge)
						{
							for (int i = 0; i < 3; i++)
							{
								double alpha = 0;
								if (uncutF(oppface, i) == oppvert)
								{
									alpha = 1.0 / 8.0;
								}
								else
								{
									alpha = 7.0 / 16.0;
								}
								periodicScoeffs.push_back({ 3 * faceitr + vertitr, 3 * oppface + i, alpha });								
							}					
						}
						else
						{
							for (int k = 0; k < 2; k++)
							{
								int target = mesh.edgeVertex(edgeitr, k);
								for (int i = 0; i < 3; i++)
								{
									if (uncutF(oppface, i) == target)
									{
										periodicScoeffs.push_back({ 3 * faceitr + vertitr, 3 * oppface + i, 1.0 });
									}
								}
							}
						}
					}
					else
					{
						int vertid = uncutF(oldfaceid, corner);
						int refface = nonProblemCopy[vertid].first;
						int refvertidx = nonProblemCopy[vertid].second;
						std::map<int, int> soupjumps;
						neighborhoodSoupJumps(mesh, problemFaces, vertexEdges, periodJumps, refface, refvertidx, soupjumps, clusterJumps);

						int eidx = 0;
						int oppface = mesh.edgeFace(edgeitr, 0);
						if (oppface == oldfaceid)
						{
							oppface = mesh.edgeFace(edgeitr, 1);
							eidx = 1;
						}						
						if (oppface == -1 || problemFaces.count(oppface))
						{
							periodicScoeffs.push_back({ 3 * faceitr + vertitr, 3 * refface + refvertidx, 1.0});
						}
						else
						{
							int opv = mesh.edgeOppositeVertex(edgeitr, eidx);
							for (int i = 0; i < 3; i++)
							{
								double alpha = 0;
								if (uncutF(oppface, i) == opv)
								{
									alpha = 1.0 / 8.0;
								}
								else
								{
									alpha = 7.0 / 16.0;
								}
								periodicScoeffs.push_back({ 3 * faceitr + vertitr, 3 * oppface + i, alpha });
								periodicb[3 * faceitr + vertitr] += -double(soupjumps[oppface]) * 2.0 * PI * alpha;								
							}							
						}
					}
				}
				else
				{
					int face1 = mesh.edgeFace(edgeitr, 0);
					int face2 = mesh.edgeFace(edgeitr, 1);
					int vert1 = mesh.edgeVertex(edgeitr, 0);
					int vert2 = mesh.edgeVertex(edgeitr, 1);

					if (face1 == -1 || face2 == -1)
					{
						int okface = face1;
						if (okface == -1)
							okface = face2;
						for (int i = 0; i < 3; i++)
						{
							if (uncutF(okface, i) == vert1 || uncutF(okface, i) == vert2)
							{
								periodicScoeffs.push_back({ 3 * faceitr + vertitr, 3 * okface + i, 1.0 / 2.0 });
							}
						}
					}
					else
					{
						int face1jmp = 0;
						int face2jmp = 0;
						if (face1 == oldfaceid)
						{
							face2jmp = periodJumps[edgeitr];
						}
						else
						{
							face1jmp = -periodJumps[edgeitr];
						}

						int oppvert1 = mesh.edgeOppositeVertex(edgeitr, 0);
						int oppvert2 = mesh.edgeOppositeVertex(edgeitr, 1);
						double beta = 1.0 / 8.0;
						int oppvalence = 0;
						int edgevalence = 0;
						if (!problemFaces.count(face1))
						{
							oppvalence++;
							edgevalence += 2;
						}
						if (!problemFaces.count(face2))
						{
							oppvalence++;
							edgevalence += 2;
						}
						double alpha = (1.0 - oppvalence * beta) / double(edgevalence);

						for (int i = 0; i < 3; i++)
						{
							if (!problemFaces.count(face1))
							{
								if (uncutF(face1, i) == oppvert1)
								{
									periodicScoeffs.push_back({ 3 * faceitr + vertitr, 3 * face1 + i, beta });
									periodicb[3 * faceitr + vertitr] += -double(face1jmp) * 2.0 * PI * beta;
								}
								else
								{
									periodicScoeffs.push_back({ 3 * faceitr + vertitr, 3 * face1 + i, alpha });
									periodicb[3 * faceitr + vertitr] += -double(face1jmp) * 2.0 * PI * alpha;
								}
							}
							if (!problemFaces.count(face2))
							{
								if (uncutF(face2, i) == oppvert2)
								{
									periodicScoeffs.push_back({ 3 * faceitr + vertitr, 3 * face2 + i, beta });
									periodicb[3 * faceitr + vertitr] += -double(face2jmp) * 2.0 * PI * beta;
								}
								else
								{
									periodicScoeffs.push_back({ 3 * faceitr + vertitr, 3 * face2 + i, alpha });
									periodicb[3 * faceitr + vertitr] += -double(face2jmp) * 2.0 * PI * alpha;
								}
							}
						}
					}
				}
			}
		}
	}

	periodicS.resize(3 * newfaces, 3 * nfaces);	
	periodicS.setFromTriplets(periodicScoeffs.begin(), periodicScoeffs.end());
}

static void zeroProblemAmplitudes(const Eigen::MatrixXi& F, const std::set<int>& problemFaces, Eigen::VectorXd& soupAmps)
{
	std::set<int> zeroverts;
	int nfaces = F.rows();
	for (int i = 0; i < nfaces; i++)
	{
		if (problemFaces.count(i) == 0)
			continue;
		for (int j = 0; j < 3; j++)
			zeroverts.insert(F(i, j));
	}

	for (int i = 0; i < nfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (zeroverts.count(F(i, j)))
			{
				soupAmps[3 * i + j] = 0;
			}
		}
	}
}

void findSharpCorners(const Eigen::MatrixXd& uncutV, const Eigen::MatrixXi& uncutF, std::set<int>& cornerVerts)
{
	int nverts = uncutV.rows();
	std::vector<double> anglesum(nverts);
	std::set<int> bdryverts;
	MeshConnectivity mesh(uncutF);
	for (int i = 0; i < uncutF.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int edge = mesh.faceEdge(i, j);
			int orient = mesh.faceEdgeOrientation(i, j);
			int opp = mesh.edgeFace(edge, 1 - orient);
			if (opp == -1)
			{
				bdryverts.insert(mesh.edgeVertex(edge, 0));
				bdryverts.insert(mesh.edgeVertex(edge, 1));
			}

			Eigen::Vector3d v0 = uncutV.row(uncutF(i, j)).transpose();
			Eigen::Vector3d v1 = uncutV.row(uncutF(i, (j + 1) % 3)).transpose();
			Eigen::Vector3d v2 = uncutV.row(uncutF(i, (j + 2) % 3)).transpose();

			Eigen::Vector3d e1 = v1 - v0;
			Eigen::Vector3d e2 = v2 - v0;
			double theta = std::atan2(e1.cross(e2).norm(), e1.dot(e2));
			anglesum[uncutF(i, j)] += theta;
		}
	}

	cornerVerts.clear();
	for (auto i : bdryverts)
	{		
		double tol = PI/4;
		if (std::fabs(anglesum[i] - PI) > tol)
			cornerVerts.insert(i);
	}
}

void wrinkledMeshUpsamplingUncut(const Eigen::MatrixXd &uncutV, const Eigen::MatrixXi &uncutF, 
	const Eigen::MatrixXd &restV, const Eigen::MatrixXi &restF,
	const Eigen::MatrixXd &cutV, const Eigen::MatrixXi &cutF, 
	const std::set<int> &noPhiFaces,
	const std::set<int> &clampedVerts,
	Eigen::MatrixXd *wrinkledV, Eigen::MatrixXi *wrinkledF, 
	Eigen::MatrixXd *upsampledTFTV, Eigen::MatrixXi *upsampledTFTF,
	Eigen::MatrixXd *soupPhiV, Eigen::MatrixXi *soupPhiF,
	Eigen::MatrixXd *soupProblemV, Eigen::MatrixXi *soupProblemF,
	Eigen::VectorXd *upsampledAmp, Eigen::VectorXd *soupPhi,
	const Eigen::VectorXd &cutAmplitude, const Eigen::VectorXd &cutPhi, 
	const SecondFundamentalFormDiscretization &sff,
	double YoungsModulus, double PoissonRatio,
	int numSubdivs, SubdivisionType subType,
	bool isUseV1Term, bool isUseV2Term)
{	
	// fix the input amplitudes
	Eigen::VectorXd fixedCutAmplitude;
	fixedCutAmplitude = cutAmplitude;

	double lameAlpha = YoungsModulus * PoissonRatio / (1.0 - PoissonRatio * PoissonRatio);
	double lameBeta = YoungsModulus / 2.0 / (1.0 + PoissonRatio);

	// find corner vertices
	std::set<int> cornerVerts;
	if (clampedVerts.size() == 0)
	{
		findSharpCorners(uncutV, uncutF, cornerVerts);
		std::cout << "Found " << cornerVerts.size() << " sharp corners" << std::endl;
	}
	else
	{
		cornerVerts = clampedVerts;
	}
	

	// turn amp and phi into a triangle soup
	Eigen::VectorXd uncutsoupamp(3 * uncutF.rows());
	Eigen::VectorXd uncutsoupphi(3 * uncutF.rows());
	for (int i = 0; i < uncutF.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			uncutsoupamp[3 * i + j] = fixedCutAmplitude[cutF(i, j)];
			uncutsoupphi[3 * i + j] = cutPhi[cutF(i, j)];
		}
	}

	// Fix phi on problem faces
	std::vector<double> canonicalPhi(uncutV.rows());
	for (int i = 0; i < uncutF.rows(); i++)
	{
		if (!noPhiFaces.count(i))
		{
			for (int j = 0; j < 3; j++)
			{
				canonicalPhi[uncutF(i, j)] = cutPhi[cutF(i, j)];
			}
		}
	}
	for (int i = 0; i < uncutF.rows(); i++)
	{
		if (noPhiFaces.count(i))
		{
			for (int j = 0; j < 3; j++)
			{
				uncutsoupphi[3 * i + j] = canonicalPhi[uncutF(i, j)];
			}
		}
	}

	std::set<int> problemFaces;
	Eigen::VectorXd soupamp;
	Eigen::VectorXd soupphi;
	//splitIsolatedProblemFaces(uncutV, uncutF, noPhiFaces, uncutsoupamp, uncutsoupphi, NV, NF, problemFaces, soupamp, soupphi);
	Eigen::MatrixXd NV = uncutV;
	Eigen::MatrixXi NF = uncutF;

	problemFaces = noPhiFaces;
	soupamp = uncutsoupamp;
	soupphi = uncutsoupphi;

	Eigen::MatrixXd RV = restV;
	Eigen::MatrixXi RF = restF;

	std::set<int> restCornerVerts;
	findCorners(RV, RF, uncutV, uncutF, restCornerVerts);	
	for (int i = 0; i < NF.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (cornerVerts.count(NF(i, j)))
				restCornerVerts.insert(RF(i, j));
		}
	}

	zeroProblemAmplitudes(NF, problemFaces, soupamp);

	// period jumps
	Eigen::VectorXi periodJumps;
	std::map<int, ClusterInfo> clusterJumps;
	computeEdgeJumps(NV.rows(), NF, problemFaces, soupphi, periodJumps, clusterJumps);

	// Subdivide everything

	for(int i=0; i<numSubdivs; ++i)
	{
		Eigen::MatrixXi tempF = NF;		
		Eigen::SparseMatrix<double> uncutS;		
		Eigen::SparseMatrix<double> soupS;
		Eigen::SparseMatrix<double> periodicS;
		Eigen::SparseMatrix<double> restS;
		Eigen::VectorXd periodicb;
		std::set<int> tempProblemFaces = problemFaces;
		Eigen::MatrixXi tempRF = RF;
		loopWithSeams(NV.rows(), tempF, periodJumps, clusterJumps, tempProblemFaces, cornerVerts,
			uncutS, soupS, periodicS, periodicb, NF, problemFaces);
		loopUncut(RV.rows(), tempRF, restCornerVerts, restS, RF);
		
		NV = (uncutS*NV).eval();
		soupamp = (soupS * soupamp).eval();
		soupphi = (periodicS * soupphi).eval() + periodicb;

		RV = (restS * RV).eval();

		zeroProblemAmplitudes(NF, problemFaces, soupamp);
		computeEdgeJumps(NV.rows(), NF, problemFaces, soupphi, periodJumps, clusterJumps);		
	}

	Eigen::MatrixXd NN;
	igl::per_vertex_normals(NV, NF, NN);

	int uncutfaces = NF.rows();	
	int uncutverts = NV.rows();
	std::vector<std::vector<double> > phivalues(uncutverts);
	std::vector<std::vector<double> > ampvalues(uncutverts);
	std::vector<std::vector<Eigen::Vector3d> > dphivalues(uncutverts);
	std::vector<std::vector<Eigen::Vector3d> > dphiperpvalues(uncutverts);
	std::vector<std::vector<double> > x0values(uncutverts);
	std::vector<std::vector<double> > y0values(uncutverts);
	MeshConnectivity finalRestMesh(RF);
	MeshConnectivity finalUncutMesh(NF);
	MeshGeometry restGeo(RV, finalRestMesh);

	Eigen::VectorXd extraDOFs;
	sff.initializeExtraDOFs(extraDOFs, finalUncutMesh, NV);

	Eigen::MatrixXd PD1, PD2;
	Eigen::VectorXd PV1, PV2;
	igl::principal_curvature(NV, NF, PD1, PD2, PV1, PV2);


	for (int i = 0; i < uncutfaces; i++)
	{
		if (problemFaces.count(i))
			continue;
		for (int j = 0; j < 3; j++)
		{
			phivalues[NF(i, j)].push_back(soupphi[3 * i + j]);
			ampvalues[NF(i, j)].push_back(soupamp[3 * i + j]);
		}
		Eigen::Matrix2d abar = firstFundamentalForm(finalRestMesh, RV, i, NULL, NULL);
		Eigen::Matrix2d a = firstFundamentalForm(finalUncutMesh, NV, i, NULL, NULL);
		
		Eigen::Matrix2d b = sff.secondFundamentalForm(finalUncutMesh, NV, extraDOFs, i, NULL, NULL);

		
		Eigen::Matrix<double, 3, 2> dr;
		dr.col(0) = (NV.row(NF(i, 1)) - NV.row(NF(i, 0))).transpose();
		dr.col(1) = (NV.row(NF(i, 2)) - NV.row(NF(i, 0))).transpose();
		
		b.setZero();
		for (int j = 0; j < 3; j++)
		{
			int vert = NF(i, j);
			Eigen::Matrix2d D;
			D << PV1[vert], 0,
				0, PV2[vert];
			Eigen::Matrix<double,2,3> U;
			U.row(0) = PD1.row(vert);
			U.row(1) = PD2.row(vert);

			b += dr.transpose() * U.transpose() * D * U * dr;
		}
		b *= 1.0 / 3.0;
		
		double dphi0 = soupphi[3 * i + 1] - soupphi[3 * i + 0];
		double dphi1 = soupphi[3 * i + 2] - soupphi[3 * i + 0];
		Eigen::Vector2d dphi(dphi0, dphi1);
		Eigen::Vector3d extdphi = dr * a.inverse() * dphi;
		for(int j=0; j<3; j++)
			dphivalues[NF(i, j)].push_back(extdphi);

		Eigen::Vector2d u = abar.inverse() * dphi;
		Eigen::Vector2d uperp = restGeo.Js.block<2,2>(2 * i, 0) * u;

		double unormsq = u.transpose() * abar * u;
		double uTbu = u.transpose() * b * u;

		double upTbup = uperp.transpose() * b * uperp;		
		double x0num = (lameAlpha / 2.0 * (abar.inverse() * b).trace() * unormsq + lameBeta * uTbu);

		double x0denom = (lameAlpha / 2.0 + lameBeta) * unormsq * unormsq;
		double x0 = unormsq > 1e-6 ? x0num / x0denom : 0;

		double y0num = 2.0 * u.transpose() * b * uperp;
		double y0denom = unormsq * unormsq;
		double y0 = unormsq > 1e-6 ? y0num / y0denom : 0;

		for (int j = 0; j < 3; j++)
		{
			x0values[NF(i, j)].push_back(x0);
			y0values[NF(i, j)].push_back(y0);
		}

		Eigen::Vector3d extdphiperp = dr * a.inverse() * abar * uperp;
		if (unormsq < 1e-6)
			extdphiperp.setZero();
		

		for(int j=0; j<3; j++)
			dphiperpvalues[NF(i, j)].push_back(extdphiperp);	
	}


	Eigen::VectorXd finalCosPhi(uncutverts);
	Eigen::VectorXd finalSin2Phi(uncutverts);
	Eigen::VectorXd finalSinPhi(uncutverts);
	Eigen::VectorXd finalAmp(uncutverts);
	Eigen::MatrixXd finalDphi(uncutverts, 3);
	Eigen::MatrixXd finalDphiperp(uncutverts, 3);
	Eigen::VectorXd finalx0(uncutverts);
	Eigen::VectorXd finaly0(uncutverts);

	double phivariance = 0;
	double ampvariance = 0;
	for (int i = 0; i < uncutverts; i++)
	{
		double meancosphi = 0;
		double meansin2phi = 0;
		double meansinphi = 0;
		for (auto it : phivalues[i])
		{
			meancosphi += std::cos(it);
			meansin2phi += std::sin(2*it);
			meansinphi += std::sin(it);
		}
		
		if (phivalues[i].size() > 0)
		{
			meancosphi /= double(phivalues[i].size());
			meansin2phi /= double(phivalues[i].size());
			meansinphi /= double(phivalues[i].size());
		}
		for (auto it : phivalues[i])
		{
			double variance = (std::cos(it) - meancosphi) * (std::cos(it) - meancosphi);
			phivariance += variance;
		}

		Eigen::Vector3d meandphi(0, 0, 0);
		Eigen::Vector3d meandphiperp(0, 0, 0);
		for (auto it : dphivalues[i])
		{
			meandphi += it;
		}
		if (dphivalues[i].size() > 0)
			meandphi /= dphivalues[i].size();
		for (auto it : dphiperpvalues[i])
		{
			meandphiperp += it;
		}
		if (dphivalues[i].size() > 0)
			meandphiperp /= dphiperpvalues[i].size();


		double meanamp = 0;
		for (auto it : ampvalues[i])
			meanamp += it;			
		
		if (ampvalues[i].size() > 0)
			meanamp /= double(ampvalues[i].size());
		for (auto it : ampvalues[i])
			ampvariance += (it - meanamp) * (it - meanamp);

		double meanx0 = 0;
		for (auto it : x0values[i])
			meanx0 += it;			

		if (x0values[i].size() > 0)
			meanx0 /= double(x0values[i].size());

		double meany0 = 0;
		for (auto it : y0values[i])
			meany0 += it;			

		if (y0values[i].size() > 0)
			meany0 /= double(y0values[i].size());

		finalCosPhi[i] = meancosphi;
		finalSin2Phi[i] = meansin2phi;
		finalSinPhi[i] = meansinphi;
		finalAmp[i] = meanamp;
		finalDphi.row(i) = meandphi.transpose();
		finalDphiperp.row(i) = meandphiperp.transpose();
		finalx0[i] = meanx0;
		finaly0[i] = meany0;
	}
	std::cout << "Variances (should be zero): " << phivariance << " " << ampvariance << std::endl;	

	if (upsampledTFTV)
		*upsampledTFTV = NV;
	if (upsampledTFTF)
		*upsampledTFTF = NF;

	int nonzerofaces = NF.rows() - problemFaces.size();
	if (soupPhiV)
		soupPhiV->resize(3 * nonzerofaces, 3);
	if (soupPhiF)
		soupPhiF->resize(nonzerofaces, 3);
	if (soupPhi)
		soupPhi->resize(3 * nonzerofaces);

	if (soupProblemV)
		soupProblemV->resize(3 * problemFaces.size(), 3);
	if (soupProblemF)
		soupProblemF->resize(problemFaces.size(), 3);

	int idx = 0;
	int pidx = 0;
	for (int i = 0; i < NF.rows(); i++)
	{
		if (problemFaces.count(i))
		{
			if (soupProblemV)
			{
				for (int j = 0; j < 3; j++)
					soupProblemV->row(3 * pidx + j) = NV.row(NF(i, j));
			}
			if (soupProblemF)
			{
				for (int j = 0; j < 3; j++)
				{
					(*soupProblemF)(pidx, j) = 3 * pidx + j;
				}
			}
			pidx++;
		}
		else
		{
			if (soupPhiV)
			{
				for (int j = 0; j < 3; j++)
				{
					soupPhiV->row(3 * idx + j) = NV.row(NF(i, j));
				}
			}
			if (soupPhiF)
			{
				for (int j = 0; j < 3; j++)
				{
					(*soupPhiF)(idx, j) = 3 * idx + j;
				}
			}
			if (soupPhi)
			{
				for (int j = 0; j < 3; j++)
				{
					(*soupPhi)[3 * idx + j] = soupphi[3 * i + j];
				}
			}
			idx++;
		}
	}

	for (int i = 0; i < uncutverts; i++)
	{
		NV.row(i) += finalAmp[i] * finalCosPhi[i] * NN.row(i);
		if (isUseV1Term)
		{
			NV.row(i) += finalAmp[i] * finalx0[i] * finalDphi.row(i) * finalSinPhi[i];
			NV.row(i) += finalAmp[i] * finaly0[i] * finalDphiperp.row(i) * finalSinPhi[i];
		}
		if(isUseV2Term)
		{
			NV.row(i) += finalAmp[i] * finalAmp[i] / 8.0 * finalDphi.row(i) * finalSin2Phi[i];
		}
		
	}

	if (wrinkledV) 
		*wrinkledV = NV;
	if (wrinkledF)
		*wrinkledF = NF;
	if (upsampledAmp)
		*upsampledAmp = finalAmp;
}

static void buildVertexNeighboringInfo(const MeshConnectivity& mesh, int nverts, std::vector<std::vector<int>>& vertNeiEdges, std::vector<std::vector<int>>& vertNeiFaces)
{
	vertNeiEdges.resize(nverts);
	vertNeiFaces.resize(nverts);

	int nfaces = mesh.nFaces();
	int nedges = mesh.nEdges();

	for(int i = 0; i < nfaces; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			vertNeiFaces[mesh.faceVertex(i, j)].push_back(i);
		}
	}

	for(int i = 0; i < nedges; i++)
	{
		for(int j = 0; j < 2; j++)
		{
			vertNeiEdges[mesh.edgeVertex(i, j)].push_back(i);
		}
	}
}

void getWrinkledMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const std::vector<std::complex<double>>& zvals, const std::vector<std::vector<int>>& vertNeiFaces, Eigen::MatrixXd& wrinkledV, bool isUseV1Term, bool isUseV2Term)
{
	int nverts = V.rows();
	int nfaces = F.rows();

	wrinkledV = V;
	Eigen::MatrixXd VN;
	igl::per_vertex_normals(V, F, VN);

	for (int vid = 0; vid < nverts; vid++)
	{
		wrinkledV.row(vid) += zvals[vid].real() * VN.row(vid);
	}

	for (int i = 0; i < nfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int vid = F(i, j);
			Eigen::Matrix2d Ib;
			Eigen::Matrix<double, 3, 2> drb;
			drb.col(0) = (V.row(F(i, (j + 1) % 3)) - V.row(F(i, j))).transpose();
			drb.col(1) = (V.row(F(i, (j + 2) % 3)) - V.row(F(i, j))).transpose();

			Ib = drb.transpose() * drb;

			std::complex<double> dz0 = zvals[F(i, (j + 1) % 3)] - zvals[F(i, j)];
			std::complex<double> dz1 = zvals[F(i, (j + 2) % 3)] - zvals[F(i, j)];

			Eigen::Vector2d aSqdtheta;
			aSqdtheta << (std::conj(zvals[vid]) * dz0).imag(), (std::conj(zvals[vid]) * dz1).imag();

			Eigen::Vector3d extASqdtheta = drb * Ib.inverse() * aSqdtheta;

			double theta = std::arg(zvals[vid]);

			if (isUseV2Term)    // now we only support v2 correction
			{
				wrinkledV.row(vid) += 1 / vertNeiFaces[vid].size() * (1. / 8 * std::sin(2 * theta) * extASqdtheta);
			}


		}
	}
}

void wrinkledMeshUpsampling(const Eigen::MatrixXd &baseV, const Eigen::MatrixXi &baseF,
                            const Eigen::MatrixXd &restV, const Eigen::MatrixXi &restF,
                            const Eigen::VectorXd &amp, const Eigen::VectorXd& omega, Eigen::VectorXd &phi,
                            Eigen::MatrixXd *wrinkledV, Eigen::MatrixXi *wrinkledF,
                            Eigen::MatrixXd *upsampledTFTV, Eigen::MatrixXi *upsampledTFTF,
                            Eigen::VectorXd *upsampledAmp, Eigen::VectorXd *upsampledPhi,
                            const SecondFundamentalFormDiscretization &sff,
                            double YoungsModulus, double PoissonRatio,
                            int numSubdivs, bool isFixBnd,
                            bool isUseV1Term, bool isUseV2Term)
{
	MeshConnectivity baseMesh(baseF);
	Eigen::VectorXd vertArea = getVertArea(baseV, baseMesh);
	Eigen::VectorXd faceArea = getFaceArea(baseV, baseMesh);
	Eigen::VectorXd edgeArea(baseMesh.nEdges());

	for (int i = 0; i < baseMesh.nEdges(); i++)
	{
		int f0 = baseMesh.edgeFace(i, 0);
		int f1 = baseMesh.edgeFace(i, 1);

		bool valid0 = f0 != -1;
		bool valid1 = f1 != -1;

		if (valid0)	// check all amp is zero
		{
			if (amp(baseF(f0, 0)) == 0 && amp(baseF(f0, 1)) == 0 && amp(baseF(f0, 2)) == 0)
			{
				valid0 = false;
			}
		}
		if (valid1)	// check all amp is zero
		{
			if (amp(baseF(f1, 0)) == 0 && amp(baseF(f1, 1)) == 0 && amp(baseF(f1, 2)) == 0)
			{
				valid1 = false;
			}
		}


		if (valid0 && valid1)
			edgeArea(i) = (faceArea(f0) + faceArea(f1)) / 2.;
		else if (valid0 && !valid1)
			edgeArea(i) = faceArea(f0) / 2.;
		else if (valid1 && !valid0)
			edgeArea(i) = faceArea(f1) / 2.;
		else
			edgeArea(i) = 0;
	}

	int nverts = baseV.rows();
	std::vector<std::complex<double>> zvals;
	roundZvalsFromEdgeOmega(baseMesh, omega, edgeArea, vertArea, nverts, zvals);

	phi.resize(nverts);
	for(int i = 0; i < nverts; i++)
	{
		phi[i] = std::arg(zvals[i]);
		zvals[i] = amp[i] * std::complex<double>(std::cos(phi[i]), std::sin(phi[i]));
	}

	Eigen::VectorXd edgeVec = swapEdgeVec(baseF, omega, 0);
	Mesh secMesh = convert2SecMesh(baseV, baseF);

	std::shared_ptr<ComplexLoop> complexLoopOpt = std::make_shared<ComplexLoopZuenko>();
	complexLoopOpt->setBndFixFlag(isFixBnd);
	complexLoopOpt->SetMesh(secMesh);

	Eigen::VectorXd subEdgeVec;
	std::vector<std::complex<double>> upZvals;

	complexLoopOpt->Subdivide(edgeVec, zvals, subEdgeVec, upZvals, numSubdivs);
	Mesh tmpMesh = complexLoopOpt->GetMesh();

	Eigen::MatrixXd upBaseV, finalWrinkleV;
	Eigen::MatrixXi upBaseF;

	parseSecMesh(tmpMesh, upBaseV, upBaseF);

	std::vector<std::vector<int>> vertNeiEdges, vertNeiFaces;
	buildVertexNeighboringInfo(MeshConnectivity(upBaseF), upBaseV.rows(), vertNeiEdges, vertNeiFaces);

	getWrinkledMesh(upBaseV, upBaseF, upZvals, vertNeiFaces, finalWrinkleV, false, isUseV2Term);

	if(wrinkledV)
	{
		*wrinkledV = finalWrinkleV;
	}
	if(wrinkledF)
	{
		*wrinkledF = upBaseF;
	}
	if(upsampledTFTV)
	{
		*upsampledTFTV = upBaseV;
	}
	if(upsampledTFTF)
	{
		*upsampledTFTF = upBaseF;
	}

	if(upsampledAmp)
	{
		upsampledAmp->resize(upZvals.size());
		for(int i = 0; i < upZvals.size(); i++)
		{
			(*upsampledAmp)[i] = std::abs(upZvals[i]);
		}
	}

	if(upsampledPhi)
	{
		upsampledPhi->resize(upZvals.size());
		for(int i = 0; i < upZvals.size(); i++)
		{
			(*upsampledPhi)[i] = std::arg(upZvals[i]);
		}
	}

}