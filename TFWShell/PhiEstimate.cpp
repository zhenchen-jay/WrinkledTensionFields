#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <deque>
#include <queue>

#include <igl/remove_unreferenced.h>
#include <igl/writeOBJ.h>
#include <igl/adjacency_list.h>
#include <igl/cotmatrix.h>
#include <igl/boundary_loop.h>

#include <Eigen/Eigenvalues> 
#include <Eigen/SPQRSupport>
#include <Eigen/CholmodSupport>

#include "PhiEstimate.h"

#include "../MeshLib/MeshConnectivity.h"
#include "../MeshLib/MeshGeometry.h"
#include "../MeshLib/IntrinsicGeometry.h"

#include "../MIPWrapper//CoMISoWrapper.h"
#include "../CommonFunctions.h"
#include "../EigenNASOQ.h"

void findCuts(const Eigen::MatrixXi &F, std::vector<std::vector<int> > &cuts)
{
	cuts.clear();

	int nfaces = F.rows();

	if (nfaces == 0)
		return;

	std::map<std::pair<int, int>, std::vector<int> > edges;
	// build edges

	for (int i = 0; i < nfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int v0 = F(i, j);
			int v1 = F(i, (j + 1) % 3);
			std::pair<int, int> e;
			e.first = std::min(v0, v1);
			e.second = std::max(v0, v1);
			edges[e].push_back(i);
		}
	}

	int nedges = edges.size();
	Eigen::MatrixXi edgeVerts(nedges,2);
	Eigen::MatrixXi edgeFaces(nedges,2);
	Eigen::MatrixXi faceEdges(nfaces, 3);
	std::set<int> boundaryEdges;
	std::map<std::pair<int, int>, int> edgeidx;
	int idx = 0;
	for (auto it : edges)
	{
		edgeidx[it.first] = idx;
		edgeVerts(idx, 0) = it.first.first;
		edgeVerts(idx, 1) = it.first.second;
		edgeFaces(idx, 0) = it.second[0];
		if (it.second.size() > 1)
		{
			edgeFaces(idx, 1) = it.second[1];
		}
		else
		{
			edgeFaces(idx, 1) = -1;
			boundaryEdges.insert(idx);
		}
		idx++;
	}
	for (int i = 0; i < nfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int v0 = F(i, j);
			int v1 = F(i, (j + 1) % 3);
			std::pair<int, int> e;
			e.first = std::min(v0, v1);
			e.second = std::max(v0, v1);
			faceEdges(i, j) = edgeidx[e];
		}
	}

	bool *deleted = new bool[nfaces];
	for (int i = 0; i < nfaces; i++)
		deleted[i] = false;

	std::set<int> deletededges;

	// loop over faces
	for (int face = 0; face < nfaces; face++)
	{
		// stop at first undeleted face
		if (deleted[face])
			continue;
		deleted[face] = true;
		std::deque<int> processEdges;
		for (int i = 0; i < 3; i++)
		{
			int e = faceEdges(face, i);
			if (boundaryEdges.count(e))
				continue;
			int ndeleted = 0;
			if (deleted[edgeFaces(e, 0)])
				ndeleted++;
			if (deleted[edgeFaces(e, 1)])
				ndeleted++;
			if (ndeleted == 1)
				processEdges.push_back(e);
		}
		// delete all faces adjacent to edges with exactly one adjacent face
		while (!processEdges.empty())
		{
			int nexte = processEdges.front();
			processEdges.pop_front();
			int todelete = -1;
			if (!deleted[edgeFaces(nexte, 0)])
				todelete = edgeFaces(nexte, 0);
			if (!deleted[edgeFaces(nexte, 1)])
				todelete = edgeFaces(nexte, 1);
			if (todelete != -1)
			{
				deletededges.insert(nexte);
				deleted[todelete] = true;
				for (int i = 0; i < 3; i++)
				{
					int e = faceEdges(todelete, i);
					if (boundaryEdges.count(e))
						continue;
					int ndeleted = 0;
					if (deleted[edgeFaces(e, 0)])
						ndeleted++;
					if (deleted[edgeFaces(e, 1)])
						ndeleted++;
					if (ndeleted == 1)
						processEdges.push_back(e);
				}
			}
		}
	}
	delete[] deleted;

	// accumulated non-deleted edges
	std::vector<int> leftedges;
	for (int i = 0; i < nedges; i++)
	{
		if (!deletededges.count(i))
			leftedges.push_back(i);
	}

	deletededges.clear();
	// prune spines
	std::map<int, std::vector<int> > spinevertedges;
	for (int i : leftedges)
	{
		spinevertedges[edgeVerts(i, 0)].push_back(i);
		spinevertedges[edgeVerts(i, 1)].push_back(i);
	}

	std::deque<int> vertsProcess;
	std::map<int, int> spinevertnbs;
	for (auto it : spinevertedges)
	{
		spinevertnbs[it.first] = it.second.size();
		if (it.second.size() == 1)
			vertsProcess.push_back(it.first);
	}
	while (!vertsProcess.empty())
	{
		int vert = vertsProcess.front();
		vertsProcess.pop_front();
		for (int e : spinevertedges[vert])
		{
			if (!deletededges.count(e))
			{
				deletededges.insert(e);
				for (int j = 0; j < 2; j++)
				{
					spinevertnbs[edgeVerts(e, j)]--;
					if (spinevertnbs[edgeVerts(e, j)] == 1)
					{
						vertsProcess.push_back(edgeVerts(e, j));
					}
				}
			}
		}
	}
	std::vector<int> loopedges;
	for (int i : leftedges)
		if (!deletededges.count(i))
			loopedges.push_back(i);

	int nloopedges = loopedges.size();
	if (nloopedges == 0)
		return;

	std::map<int, std::vector<int> > loopvertedges;
	for (int e : loopedges)
	{
		loopvertedges[edgeVerts(e, 0)].push_back(e);
		loopvertedges[edgeVerts(e, 1)].push_back(e);
	}

	std::set<int> usededges;
	for (int e : loopedges)
	{
		// make a cycle or chain starting from this edge
		while (!usededges.count(e))
		{
			std::vector<int> cycleverts;
			std::vector<int> cycleedges;
			cycleverts.push_back(edgeVerts(e, 0));
			cycleverts.push_back(edgeVerts(e, 1));
			cycleedges.push_back(e);

			std::map<int, int> cycleidx;
			cycleidx[cycleverts[0]] = 0;
			cycleidx[cycleverts[1]] = 1;

			int curvert = edgeVerts(e, 1);
			int cure = e;
			bool foundcycle = false;
			while (curvert != -1 && !foundcycle)
			{
				int nextvert = -1;
				int nexte = -1;
				for (int cande : loopvertedges[curvert])
				{
					if (!usededges.count(cande) && cande != cure)
					{
						int vidx = 0;
						if (curvert == edgeVerts(cande, vidx))
							vidx = 1;
						nextvert = edgeVerts(cande, vidx);
						nexte = cande;
						break;
					}
				}
				if (nextvert != -1)
				{
					auto it = cycleidx.find(nextvert);
					if (it != cycleidx.end())
					{
						// we've hit outselves
						std::vector<int> cut;
						for (int i = it->second; i < cycleverts.size(); i++)
						{
							cut.push_back(cycleverts[i]);
						}
						cut.push_back(nextvert);
						cuts.push_back(cut);
						for (int i = it->second; i < cycleedges.size(); i++)
						{
							usededges.insert(cycleedges[i]);
						}
						usededges.insert(nexte);
						foundcycle = true;
					}
					else
					{
						cycleidx[nextvert] = cycleverts.size();
						cycleverts.push_back(nextvert);
						cycleedges.push_back(nexte);                        
					}
				}                
				curvert = nextvert;
				cure = nexte;
			}
			if (!foundcycle)
			{
				// we've hit a dead end. reverse and try the other direction
				std::reverse(cycleverts.begin(), cycleverts.end());
				std::reverse(cycleedges.begin(), cycleedges.end());
				cycleidx.clear();
				for (int i = 0; i < cycleverts.size(); i++)
				{
					cycleidx[cycleverts[i]] = i;
				}
				curvert = cycleverts.back();
				cure = cycleedges.back();
				while (curvert != -1 && !foundcycle)
				{
					int nextvert = -1;
					int nexte = -1;
					for (int cande : loopvertedges[curvert])
					{
						if (!usededges.count(cande) && cande != cure)
						{
							int vidx = 0;
							if (curvert == edgeVerts(cande, vidx))
								vidx = 1;
							nextvert = edgeVerts(cande, vidx);
							nexte = cande;
							break;
						}
					}
					if (nextvert != -1)
					{
						auto it = cycleidx.find(nextvert);
						if (it != cycleidx.end())
						{
							// we've hit outselves
							std::vector<int> cut;
							for (int i = it->second; i < cycleverts.size(); i++)
							{
								cut.push_back(cycleverts[i]);
							}
							cut.push_back(nextvert);
							cuts.push_back(cut);
							for (int i = it->second; i < cycleedges.size(); i++)
							{
								usededges.insert(cycleedges[i]);
							}
							usededges.insert(nexte);
							foundcycle = true;
						}
						else
						{
							cycleidx[nextvert] = cycleverts.size();
							cycleverts.push_back(nextvert);
							cycleedges.push_back(nexte);                        
						}
					}                
					curvert = nextvert;
					cure = nexte;
				}
				if (!foundcycle)
				{
					// we've found a chain
					std::vector<int> cut;
					for (int i = 0; i < cycleverts.size(); i++)
					{
						cut.push_back(cycleverts[i]);
					}
					cuts.push_back(cut);
					for (int i = 0; i < cycleedges.size(); i++)
					{
						usededges.insert(cycleedges[i]);
					}                    
				}
			}
		}
	}
}


void cutMesh(const Eigen::MatrixXi &F,
	// list of cuts, each of which is a list (in order) of vertex indices of one cut.
	// Cuts can be closed loops (in which case the last vertex index should equal the
	// first) or open (in which case the two endpoint vertices should be distinct).
	// Multiple cuts can cross but there may be strange behavior if cuts share endpoint
	// vertices, or are non-edge-disjoint.
	const std::vector<std::vector<int> > &cuts,
	// new vertices and faces
	// **DO NOT ALIAS V OR F!**
	Eigen::MatrixXi &newF
)
{
	int ncuts = (int)cuts.size();

	// junction vertices that lie on multiple cuts
	std::set<int> junctions;
	std::set<int> seenverts;
	for (int i = 0; i < ncuts; i++)
	{
		std::set<int> seenincut;
		for (int j = 0; j < cuts[i].size(); j++)
		{
			if (seenverts.count(cuts[i][j]))
				junctions.insert(cuts[i][j]);
			seenincut.insert(cuts[i][j]);
		}
		for (int v : seenincut)
			seenverts.insert(v);
	}

	// "interior" cut vertices: vertices that are part of a cut but not a cut endpoint
	// or junction vertex
	std::vector<std::set<int> > cutints;
	cutints.resize(ncuts);
	for (int i = 0; i < ncuts; i++)
	{
		if (cuts[i].empty())
			continue;
		if (cuts[i].front() == cuts[i].back())
		{
			// closed loop
			for (int v : cuts[i])
			{
				if(!junctions.count(v))
					cutints[i].insert(v);                
			}
		}
		else
		{
			// open cut
			for (int j = 1; j < cuts[i].size() - 1; j++)
			{
				if(!junctions.count(cuts[i][j]))
					cutints[i].insert(cuts[i][j]);                
			}
		}
	}

	struct edge
	{
		std::pair<int, int> verts;
		edge(int v1, int v2)
		{
			verts.first = std::min(v1, v2);
			verts.second = std::max(v1, v2);
		}

		bool operator<(const edge &other) const
		{
			return verts < other.verts;
		}
	};

	// maps each edge to incident triangles
	std::map<edge, std::vector<int> > edgeTriangles;
	for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			edge e(F(i, j), F(i, (j + 1) % 3));
			edgeTriangles[e].push_back(i);
		}
	}

	// have we visited this face yet?
	bool *visited = new bool[F.rows()];

	// edges that form part of a cut
	std::set<edge> forbidden;
	for (int i = 0; i < ncuts; i++)
	{        
		for (int j = 0; j < (int)cuts[i].size() - 1; j++)
		{
			// works for both open and closed curves
			edge e(cuts[i][j], cuts[i][j + 1]);
			forbidden.insert(e);
		}
	}

	// connected components of faces adjacent to the cuts
	std::vector<std::vector<std::vector<int> > > components;
	components.resize(ncuts);

	// for each cut
	for (int cut = 0; cut < ncuts; cut++)
	{
		for (int i = 0; i < (int)F.rows(); i++)
			visited[i] = false;

		// find a face we haven't visited yet
		for (int i = 0; i < F.rows(); i++)
		{
			if (visited[i]) continue;
			bool found = false;
			for (int j = 0; j < 3; j++)
			{
				if (cutints[cut].count(F(i, j)))
				{
					found = true;
				}
			}

			if (found)
			{
				// run a BFS along the cut edges, accumulating one connected component
				// cross only edges that contain a vertex in cutints[cut], but are not forbidden
				std::deque<int> q;
				std::vector<int> component;
				q.push_back(i);
				while (!q.empty())
				{
					int next = q.front();
					q.pop_front();
					if (visited[next])
						continue;
					visited[next] = true;
					component.push_back(next);
					for (int j = 0; j < 3; j++)
					{
						int v1 = F(next, j);
						int v2 = F(next, (j + 1) % 3);
						edge e(v1, v2);
						if (cutints[cut].count(v1) == 0 && cutints[cut].count(v2) == 0)
						{
							continue;
						}
						if (forbidden.count(e))
							continue;
						for (int nb : edgeTriangles[e])
						{
							if (!visited[nb])
							{
								q.push_back(nb);
							}
						}
					} // end BFS
				}
				components[cut].push_back(component);
			} // end if found
		} // end loop over all faces
	} // end loop over cuts

	std::map<int, std::vector<std::vector<int> > > junctcomponents;

	// for each junction
	for (int junc : junctions)
	{
		for (int i = 0; i < (int)F.rows(); i++)
			visited[i] = false;

		// find a face we haven't visited yet
		for (int i = 0; i < F.rows(); i++)
		{
			if (visited[i]) continue;
			bool found = false;
			for (int j = 0; j < 3; j++)
			{
				if (junc == F(i, j))
				{
					found = true;
				}
			}

			if (found)
			{
				// run a BFS along the cut edges, accumulating one connected component
				// cross only edges that contain the junction, but are not forbidden
				std::deque<int> q;
				std::vector<int> component;
				q.push_back(i);
				while (!q.empty())
				{
					int next = q.front();
					q.pop_front();
					if (visited[next])
						continue;
					visited[next] = true;
					component.push_back(next);
					for (int j = 0; j < 3; j++)
					{
						int v1 = F(next, j);
						int v2 = F(next, (j + 1) % 3);
						edge e(v1, v2);
						if (v1 != junc && v2 != junc)
						{
							continue;
						}
						if (forbidden.count(e))
							continue;
						for (int nb : edgeTriangles[e])
						{
							if (!visited[nb])
							{
								q.push_back(nb);
							}
						}
					} // end BFS
				}
				junctcomponents[junc].push_back(component);
			} // end if found
		} // end loop over all faces
	} // end loop over cuts

	int vertstoadd = 0;
	// create a copy of each vertex for each component of each cut
	for (int i = 0; i < ncuts; i++)
	{
		vertstoadd += components[i].size() * cutints[i].size();
	}
	// create a copy of each junction point for each component of each junction
	for (int v : junctions)
	{
		vertstoadd += junctcomponents[v].size();
	}

	// create new faces
	Eigen::MatrixXi augF = F;

	// duplicate vertices and reindex faces

	int idx = 0;
	for (int i = 0; i < F.rows(); i++)
		for (int j = 0; j < 3; j++)
			idx = std::max(idx, F(i, j));
	idx++;

	for (int cut = 0; cut < ncuts; cut++)
	{
		for (int i = 0; i < components[cut].size(); i++)
		{
			// duplicate vertices
			std::map<int, int> idxmap;
			for (int v : cutints[cut])
			{
				idxmap[v] = idx;
				idx++;
			}
			for (int f : components[cut][i])
			{
				for (int j = 0; j < 3; j++)
				{
					int v = augF(f, j);
					if (cutints[cut].count(v))
						augF(f, j) = idxmap[v];
				}
			}
		}        
	}

	for (int junc : junctions)
	{
		for (int i = 0; i < junctcomponents[junc].size(); i++)
		{

			int newidx = idx;
			idx++;

			for (int f : junctcomponents[junc][i])
			{
				for (int j = 0; j < 3; j++)
				{
					int v = augF(f, j);
					if (v == junc)
						augF(f, j) = newidx;
				}
			}
		}
	}

	newF = augF;
	delete[] visited;
}

static double angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, const Eigen::Vector3d axis)
{
	return 2.0 * atan2(v1.cross(v2).dot(axis), v1.norm() * v2.norm() + v1.dot(v2));
}

void vectorFieldSingularities(const Eigen::MatrixXi &F, const std::vector<Eigen::Matrix2d> &abars, const Eigen::MatrixXd &w, std::vector<int> &singularities)
{
	singularities.clear();
	MeshConnectivity mesh(F);
	IntrinsicGeometry geom(mesh, abars);
	
	std::set<int> checked;

	int nfaces = F.rows();
	for(int i=0; i<nfaces; i++)
	{
		for(int j=0; j<3; j++)
		{
			int centervert = F(i,j);
			if(!checked.count(centervert))
			{
				int startface = i;
				int startspoke = (j+1)%3;
		
				double wangle = 0;

				int curface = startface;
				int curspoke = startspoke;
				double totangle = 0;

				bool isboundary = false;

				while (true)
				{
					int edge = mesh.faceEdge(curface, curspoke);
					int side = (mesh.edgeFace(edge, 0) == curface) ? 0 : 1;
					int nextface = mesh.edgeFace(edge, 1 - side);
					if (nextface == -1)
					{
						isboundary = true;
						break;
					}
			
					Eigen::Vector2d curw = abars[curface].inverse() * w.row(curface).transpose();
					Eigen::Vector2d nextwbary = abars[nextface].inverse() * w.row(nextface).transpose();
					
					Eigen::Vector2d nextw = geom.Ts.block<2, 2>(2 * edge, 2 - 2 * side) * nextwbary;
					Eigen::Vector2d nextwperp = geom.Js.block<2, 2>(2 * curface, 0) * nextw;
					double curwnorm = curw.transpose() * abars[curface] * curw;
					double nextwnorm = nextw.transpose() * abars[curface] * nextw;
					double crossprod = (nextwperp.transpose() * abars[curface] * curw);
					double innerprod = (nextw.transpose() * abars[curface] * curw);
					double angleopt = std::atan2(crossprod, innerprod);
					
					wangle += angleopt;
					
					int spokep1 = (curspoke + 1) % 3;
					int apex = (curspoke + 2) % 3;
					Eigen::Vector2d barys[3] = { {0,0}, {1,0}, {0,1} };
					Eigen::Vector2d edge1 = barys[curspoke] - barys[apex];
					Eigen::Vector2d edge2 = barys[spokep1] - barys[apex];
					double e1norm = std::sqrt(edge1.transpose() * abars[curface] * edge1);
					double e2norm = std::sqrt(edge2.transpose() * abars[curface] * edge2);
					double eprod = edge1.transpose() * abars[curface] * edge2;
					double cose = std::min(std::max(eprod / e1norm / e2norm, -1.0), 1.0);
					totangle += std::acos(cose);

					curface = nextface;
					for (int k = 0; k < 3; k++)
					{
						if (F(nextface, k) == centervert)
						{
							curspoke = (k + 1) % 3;
							break;
						}
					}
		
					if (curface == startface)
						break;
				}

				if (!isboundary)
				{
					const double PI = 3.1415926535898;
					double index = wangle + 2 * PI - totangle;                    
					if (fabs(index) > PI)
					{
						std::cout << wangle << " " << totangle << std::endl;
						singularities.push_back(centervert);
					}
				}
				checked.insert(centervert);
			}
		}
	}
}

void punctureMesh(const Eigen::MatrixXi &F, const std::vector<int> &singularities, Eigen::MatrixXi &puncturedF, Eigen::VectorXi &newFacesToOld)
{
	std::vector<int> okfaces;
	int nfaces = F.rows();
	std::set<int> singularset;
	for(auto it : singularities)
		singularset.insert(it);
		
	for(int i=0; i<nfaces; i++)
	{
		bool ok = true;
		for(int j=0; j<3; j++)
		{
			if(singularset.count(F(i,j)))
				ok = false;
		}
		if(ok)
			okfaces.push_back(i);
	}
	
	puncturedF.resize(okfaces.size(), 3);
	newFacesToOld.resize(okfaces.size());
	int idx=0;
	for(auto it : okfaces)
	{
		newFacesToOld[idx] = it;
		puncturedF.row(idx) = F.row(it);
		idx++;
	}
}

void punctureMeshUsingPureTension(const Eigen::MatrixXi& F, const std::set<int>& tensionFaces, Eigen::MatrixXi& puncturedF, Eigen::VectorXi& newFacesToOld)
{
	std::vector<int> okfaces;
	int nfaces = F.rows();

	for (int i = 0; i < nfaces; i++)
	{
		bool ok = true;
		if (tensionFaces.count(i))
			ok = false;
		if (ok)
			okfaces.push_back(i);
	}

	puncturedF.resize(okfaces.size(), 3);
	newFacesToOld.resize(okfaces.size());
	int idx = 0;
	for (auto it : okfaces)
	{
		newFacesToOld[idx] = it;
		puncturedF.row(idx) = F.row(it);
		idx++;
	}
}

struct UnionFind
{
	std::vector<int> parent;
	std::vector<int> sign;
	UnionFind(int items) 
	{
		parent.resize(items);
		sign.resize(items);
		for(int i=0; i<items; i++)
		{
			parent[i] = i;
			sign[i] = 1;
		}
	}
		
	std::pair<int, int> find(int i)
	{
		if(parent[i] != i)
		{
			auto newparent = find(parent[i]);
			sign[i] *= newparent.second;
			parent[i] = newparent.first;
		}
		
		return {parent[i], sign[i]};
	}
	
	void dounion(int i, int j, int usign)
	{
		auto xroot = find(i);
		auto yroot = find(j);
		if(xroot.first != yroot.first)
		{
			parent[xroot.first] = yroot.first;
			sign[xroot.first] = usign * xroot.second * yroot.second;
		}
	}
};

static const double PI = 3.1415926535898;

void combField(const Eigen::MatrixXi &F, 
	const std::vector<Eigen::Matrix2d> &abars, 
	const Eigen::VectorXd* weight,
	const Eigen::MatrixXd &w, Eigen::MatrixXd &combedW)
{
	int nfaces = F.rows();
	UnionFind uf(nfaces);
	MeshConnectivity mesh(F);
	IntrinsicGeometry geom(mesh, abars);
	struct Visit
	{
		int edge;
		int sign;
		double norm;
		bool operator<(const Visit &other) const
		{
			return norm > other.norm;
		}
	};
	
	std::priority_queue<Visit> pq;
	
	int nedges = mesh.nEdges();
	for(int i=0; i<nedges; i++)
	{
		int face1 = mesh.edgeFace(i, 0);
		int face2 = mesh.edgeFace(i, 1);
		if(face1 == -1 || face2 == -1)
			continue;
			
		Eigen::Vector2d curw = abars[face1].inverse() * w.row(face1).transpose();
		Eigen::Vector2d nextwbary1 = abars[face2].inverse() * w.row(face2).transpose();                    
		Eigen::Vector2d nextw = geom.Ts.block<2, 2>(2 * i, 2) * nextwbary1;
		int sign = ( (curw.transpose() * abars[face1] * nextw) < 0 ? -1 : 1);
		double innerp = curw.transpose() * abars[face1] * nextw;

		if (!weight)
		{
			double normcw = std::sqrt(curw.transpose() * abars[face1] * curw);
			double normnw = std::sqrt(nextw.transpose() * abars[face1] * nextw);
			double negnorm = -std::min(normcw, normnw);
			pq.push({ i, sign, negnorm });

		}
		else
		{
			double compcw = weight->coeffRef(face1);
			double compnw = weight->coeffRef(face2);
			double negw = -std::min(compcw, compnw);

			pq.push({ i, sign, negw });
		}

		
	}
	
	while(!pq.empty())
	{
		auto next = pq.top();
		pq.pop();
		uf.dounion(mesh.edgeFace(next.edge, 0), mesh.edgeFace(next.edge, 1), next.sign);
	}
	
	combedW.resize(nfaces, 2);
	for(int i=0; i<nfaces; i++)
	{
		int sign = uf.find(i).second;
		combedW.row(i) = w.row(i) * sign;
	}

}

void combFieldCutbyTension(const Eigen::MatrixXi& F,
	const std::vector<Eigen::Matrix2d>& abars,
	const std::set<int> tensionFaces,
	const Eigen::MatrixXd& w, Eigen::MatrixXd& combedW)
{
	Eigen::MatrixXi puncturedF;
	Eigen::VectorXi newFacesToOld;
	punctureMeshUsingPureTension(F, tensionFaces, puncturedF, newFacesToOld);

	MeshConnectivity puncturedMesh(puncturedF);
	int nPuncturedFaces = puncturedMesh.nFaces();
	Eigen::MatrixXd puncturedW(nPuncturedFaces, 2);
	std::vector<Eigen::Matrix2d> puncturedAbars(nPuncturedFaces);

	for (int i = 0; i < nPuncturedFaces; i++)
	{
		puncturedW.row(i) = w.row(newFacesToOld(i));
		puncturedAbars[i] = abars[newFacesToOld(i)];
	}

	Eigen::MatrixXd combedPuncturedW;

	combField(puncturedF, puncturedAbars, NULL, puncturedW, combedPuncturedW);

	combedW = w;
	for (int i = 0; i < nPuncturedFaces; i++)
	{
		combedW.row(newFacesToOld(i)) = combedPuncturedW.row(i);
	}

}


void reindex(Eigen::MatrixXi &F)
{
	std::map<int, int> old2new;
	for(int i=0; i<F.rows(); i++)
	{
		for(int j=0; j<3; j++)
		{
			int vidx = F(i,j);
			if(old2new.find(vidx) == old2new.end())
			{
				int newidx = old2new.size();
				old2new[vidx] = newidx;
			}
		}
	}
	for(int i=0; i<F.rows(); i++)
	{
		for(int j=0; j<3; j++)
		{
			int vidx = F(i,j);
			int newvidx = old2new[vidx];
			F(i,j) = newvidx;
		}
	}
}


void removeLocalSingularity(const Eigen::VectorXd& phi, const Eigen::MatrixXi& F, Eigen::VectorXd& newPhi)
{
	std::vector<std::vector<int>> adjacencyLists;
	MeshConnectivity mesh(F);

	int nverts = phi.rows();
	igl::adjacency_list(F, adjacencyLists);

	newPhi = phi;


	double absAvePhiDiff = 0;
	for (int i = 0; i < mesh.nEdges(); i++)
	{
		double phiDiff = phi(mesh.edgeVertex(i, 0)) - phi(mesh.edgeVertex(i, 1));
		absAvePhiDiff += abs(phiDiff);
	}

	absAvePhiDiff /= mesh.nEdges();
	Eigen::VectorXi singularity(nverts);
	singularity.setZero();

	for (int i = 0; i < nverts; i++)
	{
		const std::vector<int> neighbors = adjacencyLists[i];
		double centerPhi = phi(i);

		bool isSingularity = true;
		int count = 0;
		for (int j = 0; j < neighbors.size(); j++)
		{
			if (abs(phi(neighbors[j]) - centerPhi) < absAvePhiDiff)
				isSingularity = false;
			else
				count++;
		}

		//if(!isSingularity) // check 2-rings
		//{
		//	for (int j = 0; j < neighbors.size(); j++)
		//	{
		//		int centerCount = count;
		//		// check 2 neighbors
		//		int neigVert = neighbors[j];
		//		double anothercenterPhi = phi(neigVert);
		//		for (int k = 0; k < adjacencyLists[neigVert].size(); k++)
		//		{
		//			if (abs(phi(adjacencyLists[neigVert][k]) - anothercenterPhi) > absAvePhiDiff)
		//				centerCount++;
		//		}
		//		if (centerCount >= neighbors.size() + adjacencyLists[neigVert].size() - 2)
		//			isSingularity = true;
		//	}
		//}
		if (isSingularity)
			singularity(i) = 1;

	}

	int numSingularity = 0;

	for (int i = 0; i < nverts; i++)
	{
		if (singularity(i))
		{
			numSingularity++;
			const std::vector<int> neighbors = adjacencyLists[i];

			newPhi(i) = 0;
		   

			for (int j = 0; j < neighbors.size(); j++)
			{
				newPhi(i) += phi(neighbors[j]);
			   
			}

			newPhi(i) /= neighbors.size();
		}
	}
	std::cout << "number of singularities in phi: " << numSingularity << std::endl;


}


void estimateAmpOmegaFromStrain(const std::vector<Eigen::Matrix2d>& abars,
	const Eigen::MatrixXd& curPos,
	const Eigen::MatrixXi& F,
	const std::set<int>& clampedVerts,
	double amplitudeEstimate,
	Eigen::VectorXd& amp,
	Eigen::MatrixXd& w,
	std::set<int>& tensionFaces)
{
	MeshConnectivity mesh(F);
	MeshGeometry curGeo(curPos, mesh);
	assert(abars.size() == F.rows());
	int nfaces = F.rows();

	Eigen::MatrixXd wguess(nfaces, 2);
	wguess.setZero();

	std::set<int> pureTensionVerts;

	double length = curPos.row(0).maxCoeff() - curPos.row(0).minCoeff();
	double width = curPos.row(1).maxCoeff() - curPos.row(1).minCoeff();
	double height = curPos.row(2).maxCoeff() - curPos.row(2).minCoeff();

	double bboxSize = std::max(length, width);
	bboxSize = std::max(bboxSize, height);
	std::cout << "bbox: " << length << ",  " << width << ", " << height << std::endl;

	locatePotentialPureTensionFaces(abars, curPos, F, tensionFaces);
	double newAmpEstimate = std::min(amplitudeEstimate, 0.1 * bboxSize);

	Eigen::VectorXd faceAmp(nfaces);
	faceAmp.setZero();

	for (int i = 0; i < nfaces; i++)
	{
		if (tensionFaces.find(i) != tensionFaces.end())
			continue;
		Eigen::Matrix2d abar = abars[i];
		Eigen::Matrix2d a = curGeo.Bs[i].transpose() * curGeo.Bs[i];
		Eigen::Matrix2d diff = a - abar;
		Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix2d> solver(diff, abar);
		Eigen::Vector2d evec = solver.eigenvectors().col(0);

		faceAmp(i) = std::sqrt(2.0 * std::fabs(solver.eigenvalues()[0]));
		wguess.row(i) = (abar * evec).transpose();	// convert to one-form
	}

	// rescale the faceAmp such that we don't have a huge amp
	double coeff = newAmpEstimate / (faceAmp.maxCoeff() + 1e-6);;
	faceAmp = faceAmp * coeff;
	wguess = wguess / coeff;

	std::cout << "w estimated finished! " << std::endl;
	std::cout << "max face amp = " << faceAmp.maxCoeff() << std::endl;

	Eigen::MatrixXd combedW = wguess;
	combField(F, abars, &faceAmp, wguess, combedW);

	w = combedW;

	amp.resize(curPos.rows());
	amp.setZero();

	Eigen::VectorXd vdegree(curPos.rows());
	vdegree.setZero();

	for (int i = 0; i < nfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int vid = mesh.faceVertex(i, j);
			vdegree(vid) += 1.0;
			amp(vid) += faceAmp(i);
		}
	}
	for (int i = 0; i < curPos.rows(); i++)
	{
		amp(i) /= vdegree(i);
	}


	for (auto& f : tensionFaces)
	{
		for (int j = 0; j < 3; j++)
		{
			int vid = mesh.faceVertex(f, j);
			if (pureTensionVerts.find(vid) == pureTensionVerts.end())
				pureTensionVerts.insert(vid);
		}
	}

	for (auto& vid : pureTensionVerts)
		amp(vid) = 0;
	for (auto& vid : clampedVerts)
		amp(vid) = 0;

	// smoothing for the amplitude
	std::vector<std::vector<int>> vertNeis;
	igl::adjacency_list(F, vertNeis);

	Eigen::SparseMatrix<double> L(curPos.rows(), curPos.rows());
	std::vector<Eigen::Triplet<double>> Lcoeff;

	for (int i = 0; i < curPos.rows(); i++)
	{
		if (amp(i) > 0) // not pure tension vertices or clamped vertices
		{
			for (auto& neiV : vertNeis[i])
			{
				Lcoeff.push_back(Eigen::Triplet<double>(i, neiV, 1.0 / vertNeis[i].size()));
			}
		}
		else
			Lcoeff.push_back(Eigen::Triplet<double>(i, i, 1.0));
	}

	L.setFromTriplets(Lcoeff.begin(), Lcoeff.end());

	// apply twice by default
	amp = L * (L * amp);
}

void roundPhiFromOmega(
	const Eigen::MatrixXi& F,
	const std::vector<Eigen::Matrix2d>& abars,
	const Eigen::MatrixXd& w,
	Eigen::VectorXd& phi,
	Eigen::MatrixXd& roundedW)
{
	int nfaces = F.rows();

	std::vector<std::vector<int> > cuts;
	findCuts(F, cuts);
	std::cout << "Used " << cuts.size() << " cuts" << std::endl;

	std::map<std::pair<int, int>, int> cutpairs;

	// convert cuts to edge indices
	for (int i = 0; i < cuts.size(); i++)
	{
		int len = cuts[i].size();
		for (int j = 0; j < len - 1; j++)
		{
			cutpairs[std::pair<int, int>(cuts[i][j], cuts[i][j + 1])] = i;
		}
	}

	Eigen::MatrixXi newF;

	cutMesh(F, cuts, newF);
	reindex(newF);

	int newverts = 0;
	for (int i = 0; i < newF.rows(); i++)
		for (int j = 0; j < 3; j++)
			newverts = std::max(newverts, newF(i, j) + 1);

	int intdofs = 0;
	MeshConnectivity punctmesh(F);

	std::vector<Eigen::Triplet<double> > Ccoeffs;

	// integer constraints on the two sides of the cut
	int row = 0;
	for (int i = 0; i < punctmesh.nEdges(); i++)
	{
		int f1 = punctmesh.edgeFace(i, 0);
		int f2 = punctmesh.edgeFace(i, 1);
		if (f1 == -1 || f2 == -1)
			continue;
		int v1 = punctmesh.edgeVertex(i, 0);
		int v2 = punctmesh.edgeVertex(i, 1);
		double sign = 1.0;
		auto it = cutpairs.find(std::pair<int, int>(v1, v2));
		if (it == cutpairs.end())
		{
			it = cutpairs.find(std::pair<int, int>(v2, v1));
			sign = -1.0;
		}
		if (it == cutpairs.end())
			continue;

		int newv1 = -1;
		int newv2 = -1;
		int neww1 = -1;
		int neww2 = -1;
		for (int j = 0; j < 3; j++)
		{
			if (F(f1, j) == v1)
				newv1 = newF(f1, j);
			if (F(f2, j) == v1)
				newv2 = newF(f2, j);
			if (F(f1, j) == v2)
				neww1 = newF(f1, j);
			if (F(f2, j) == v2)
				neww2 = newF(f2, j);
		}

		Ccoeffs.push_back(Eigen::Triplet<double>(row, newv1, sign));
		Ccoeffs.push_back(Eigen::Triplet<double>(row, newv2, -sign));
		Ccoeffs.push_back(Eigen::Triplet<double>(row, newverts + intdofs, 1.0));
		Ccoeffs.push_back(Eigen::Triplet<double>(row + 1, neww1, sign));
		Ccoeffs.push_back(Eigen::Triplet<double>(row + 1, neww2, -sign));
		Ccoeffs.push_back(Eigen::Triplet<double>(row + 1, newverts + intdofs, 1.0));
		row += 2;
		intdofs++;
	}
	Eigen::SparseMatrix<double> C(row, newverts + intdofs + 1);
	C.setFromTriplets(Ccoeffs.begin(), Ccoeffs.end());

	std::vector<Eigen::Triplet<double> > Acoeffs;
	int newfaces = newF.rows();
	//assert(newfaces == punctF.rows());

	for (int i = 0; i < newfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Acoeffs.push_back(Eigen::Triplet<double>(3 * i + j, newF(i, j), 1.0));
		}
	}
	Eigen::SparseMatrix<double> A(3 * newfaces, newverts + intdofs);
	A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

	std::vector<Eigen::Triplet<double> > L2Minvcoeffs;
	std::vector<Eigen::Triplet<double> > Mcoeffs;
	for (int i = 0; i < newfaces; i++)
	{
		double area = 0.5 * std::sqrt(abars[i].determinant());
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < 2; k++)
			{
				L2Minvcoeffs.push_back(Eigen::Triplet<double>(2 * i + j, 2 * i + k, area * abars[i].inverse()(j, k)));
				Mcoeffs.push_back(Eigen::Triplet<double>(2 * i + j, 2 * i + k, abars[i](j, k)));
			}
	}
	Eigen::SparseMatrix<double> L2Minv(2 * newfaces, 2 * newfaces);
	L2Minv.setFromTriplets(L2Minvcoeffs.begin(), L2Minvcoeffs.end());
	Eigen::SparseMatrix<double> M(2 * newfaces, 2 * newfaces);
	M.setFromTriplets(Mcoeffs.begin(), Mcoeffs.end());

	std::vector<Eigen::Triplet<double> > Dcoeffs;
	for (int i = 0; i < newfaces; i++)
	{
		Dcoeffs.push_back(Eigen::Triplet<double>(2 * i, 3 * i + 1, 1.0));
		Dcoeffs.push_back(Eigen::Triplet<double>(2 * i, 3 * i, -1.0));
		Dcoeffs.push_back(Eigen::Triplet<double>(2 * i + 1, 3 * i + 2, 1.0));
		Dcoeffs.push_back(Eigen::Triplet<double>(2 * i + 1, 3 * i, -1.0));
	}
	Eigen::SparseMatrix<double> D(2 * newfaces, 3 * newfaces);
	D.setFromTriplets(Dcoeffs.begin(), Dcoeffs.end());

	Eigen::SparseMatrix<double> Mat = A.transpose() * D.transpose() * L2Minv * D * A;
	Eigen::VectorXd punctfield(2 * newfaces);
	for (int i = 0; i < newfaces; i++)
	{
		punctfield[2 * i] = w(i, 0) / 2.0 / PI;
		punctfield[2 * i + 1] = w(i, 1) / 2.0 / PI;
	}
	Eigen::VectorXd rhs = A.transpose() * D.transpose() * L2Minv * punctfield;

	Eigen::VectorXd result;
	Eigen::VectorXi toRound(intdofs);
	for (int i = 0; i < intdofs; i++)
	{
		toRound[i] = newverts + i;
	}
	std::cout << intdofs << " integer variables" << std::endl;
	ComisoWrapper(C, Mat, result, rhs, toRound, 1e-6);
	std::cout << "Solver residual: " << (Mat * result - rhs).norm() << std::endl;
	std::cout << "Reconstruction residual: " << (D * A * result - punctfield).transpose() * L2Minv * (D * A * result - punctfield) << std::endl;

	// map phi back to the original mesh
	int oldverts = 0;
	for (int i = 0; i < F.rows(); i++)
		for (int j = 0; j < 3; j++)
			oldverts = std::max(oldverts, F(i, j) + 1);

	phi.resize(oldverts);
	phi.setZero();
	for (int i = 0; i < newfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			phi[F(i, j)] = 2.0 * PI * result[newF(i, j)];
		}
	}

	Eigen::VectorXd cutPhi;
	cutPhi.resize(newverts);
	roundedW.resize(newfaces, 2);
	roundedW.setZero();

	for (int i = 0; i < newfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			cutPhi(newF(i, j)) = 2.0 * PI * result(newF(i, j));
		}
	}
	Eigen::VectorXd newPhi;
	removeLocalSingularity(cutPhi, newF, newPhi);

	for (int i = 0; i < newfaces; i++)
	{
		roundedW(i, 0) = newPhi(newF(i, 1)) - newPhi(newF(i, 0));
		roundedW(i, 1) = newPhi(newF(i, 2)) - newPhi(newF(i, 0));
	}

}

void estimateWrinkleVariablesFromStrain(
	const std::vector<Eigen::Matrix2d>& abars,
	const Eigen::MatrixXd& curPos,
	const Eigen::MatrixXi& F,
	const std::set<int>& clampedVerts,
	double amplitudeEstimate,
	Eigen::VectorXd& amp,
	Eigen::VectorXd& phi,
	Eigen::VectorXd& dphi,
	std::set<int>& tensionFaces)
{
	MeshConnectivity mesh(F);
	MeshGeometry curGeo(curPos, mesh);
	assert(abars.size() == F.rows());
	int nfaces = F.rows();

	Eigen::MatrixXd w;
	estimateAmpOmegaFromStrain(abars, curPos, F, clampedVerts, amplitudeEstimate, amp, w, tensionFaces);

	std::vector<int> singularVertices;
	vectorFieldSingularities(F, abars, w, singularVertices);
	std::cout << "Found " << singularVertices.size() << " singularities" << std::endl;
	Eigen::MatrixXi punctF;
	Eigen::VectorXi new2old;
	punctureMesh(F, singularVertices, punctF, new2old);
	int newfaces = punctF.rows();
	Eigen::VectorXd punctphi;

	std::vector<Eigen::Matrix2d> punctabars;
	for (int i = 0; i < newfaces; i++)
		punctabars.push_back(abars[new2old[i]]);

	Eigen::MatrixXd punctW(newfaces, 2);
	for (int i = 0; i < newfaces; i++)
		punctW.row(i) = w.row(new2old[i]);

	Eigen::MatrixXd roundedW;
	std::cout << "Round phi from omega" << std::endl;
	roundPhiFromOmega(punctF, punctabars, punctW, punctphi, roundedW);
	std::cout << "Round phi from omega finished" << std::endl;

	phi.resize(curPos.rows());
	phi.setZero();

	Eigen::MatrixXd faceDphi(F.rows(), 2);
	faceDphi.setZero();

	for (int i = 0; i < newfaces; i++)
	{
		int oldface = new2old[i];
		for (int j = 0; j < 3; j++)
		{
			phi[F(oldface, j)] = punctphi[punctF(i, j)];
		}
		faceDphi.row(oldface) = roundedW.row(i);
	}

	std::vector<std::vector<int>> adjLst;
	igl::adjacency_list(F, adjLst);
	auto backupPhi = phi;

	for (auto& it : singularVertices)
	{
		if (backupPhi(it))
		{
			std::cout << "something wrong happened" << std::endl;
		}
		for (int j = 0; j < adjLst[it].size(); j++)
		{
			phi(it) += backupPhi(adjLst[it][j]);
		}
		phi(it) /= adjLst[it].size();
	}

	Eigen::VectorXd newPhi;
	removeLocalSingularity(phi, F, newPhi);

	phi = newPhi;

	int nedges = mesh.nEdges();
	dphi.resize(nedges);
	dphi.setZero();
	for (int i = 0; i < nfaces; i++)
	{
		for (int e = 0; e < 3; e++)
		{
			int edge = mesh.faceEdge(i, e);
			int vert1 = mesh.edgeVertex(edge, 0);
			int vert2 = mesh.edgeVertex(edge, 1);
			Eigen::Vector2d barys[3] = { {0,0}, {1,0}, {0,1} };
			Eigen::Vector2d facee(0, 0);
			for (int j = 0; j < 3; j++)
			{
				if (vert1 == mesh.faceVertex(i, j))
					facee -= barys[j];
				else if (vert2 == mesh.faceVertex(i, j))
					facee += barys[j];
			}
			dphi[edge] = faceDphi.row(i).dot(facee);
		}
	}
}

void faceDPhi2EdgeDPhi(const Eigen::MatrixXd& faceDphi, 
	const std::set<int>& tensionFaces, 
	const std::vector<Eigen::Matrix2d>& abars,
	Eigen::MatrixXi F, 
	Eigen::VectorXd& dphi)
{
	MeshConnectivity mesh(F);
	int nedges = mesh.nEdges();
	int nfaces = mesh.nFaces();
	if (nfaces != faceDphi.rows())
	{
		std::cout << "Error in face dphi value" << std::endl;
	}

	Eigen::SparseMatrix<double> C;
	std::vector<Eigen::Triplet<double> > constraintCoeff;

	int row = 0;
	for (int i = 0; i < nfaces; i++)
	{
		if (tensionFaces.find(i) != tensionFaces.end())
			continue;

		for (int j = 0; j < 3; j++)
		{
			int eid = mesh.faceEdge(i, j);
			int vert1 = mesh.edgeVertex(eid, 0);
			int vert2 = mesh.edgeVertex(eid, 1);

			if (vert1 == mesh.faceVertex(i, (j + 1) % 3))
				constraintCoeff.push_back(Eigen::Triplet<double>(row, eid, 1.0));
			else
				constraintCoeff.push_back(Eigen::Triplet<double>(row, eid, -1.0));
		}

		row++;
	}
	C.resize(row, nedges);
	C.setFromTriplets(constraintCoeff.begin(), constraintCoeff.end());
	std::cout << nfaces - row << " faces are removed from integrabilty due to pure tension check." << std::endl;

	std::vector<Eigen::Triplet<double> > Mcoeffs;	
	for (int i = 0; i < nfaces; i++)
	{
		Eigen::Matrix2d ainv = abars[i].inverse();
		double area = std::sqrt(abars[i].determinant());
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				Mcoeffs.push_back({ 2 * i + j, 2 * i + k, ainv(j,k) * area });
			}
		}
	}
	Eigen::SparseMatrix<double> M(2 * nfaces, 2 * nfaces);
	M.setFromTriplets(Mcoeffs.begin(), Mcoeffs.end());

	std::vector<Eigen::Triplet<double> > L;
	Eigen::SparseMatrix<double> edgeDphi2FaceDphi(2 * nfaces, nedges);

	for (int i = 0; i < nfaces; i++)
	{
		if (tensionFaces.find(i) != tensionFaces.end())
			continue;
		Eigen::Vector3i edgeIndices;
		for (int j = 0; j < 3; j++)
		{
			edgeIndices(j) = mesh.faceEdge(i, j);
		}

		double flagU = 1.0;
		double flagV = 1.0;

		if (mesh.faceVertex(i, 0) > mesh.faceVertex(i, 1))
		{
			flagU = -1;
		}
		if (mesh.faceVertex(i, 0) > mesh.faceVertex(i, 2))
		{
			flagV = -1;
		}

		L.push_back(Eigen::Triplet<double>(2 * i, edgeIndices(2), flagU));
		L.push_back(Eigen::Triplet<double>(2 * i + 1, edgeIndices(1), flagV));
	}
	edgeDphi2FaceDphi.setFromTriplets(L.begin(), L.end());

	double reg = 1e-8;
	Eigen::VectorXd B, Beq, Bieq;
	Eigen::SparseMatrix<double> Q, Aeq, Aieq;
	Eigen::VectorXd lx, ux, delta_x;

	Q = edgeDphi2FaceDphi.transpose() * M * edgeDphi2FaceDphi;
	Aeq = C;
	Beq = Eigen::VectorXd::Zero(Aeq.rows());

	lx.resize(0);
	ux.resize(0);

	Bieq.resize(0);
	Aieq.resize(0, 0);

	Eigen::VectorXd w(2 * nfaces);
	for (int i = 0; i < nfaces; i++)
	{
		if (tensionFaces.find(i) != tensionFaces.end())
		{
			w(2 * i) = 0;
			w(2 * i + 1) = 0;
		}
		else
		{
			w(2 * i) = faceDphi(i, 0);
			w(2 * i + 1) = faceDphi(i, 1);
		}
	}
	std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << "|| w || = " << w.norm() << std::endl;
	B = -edgeDphi2FaceDphi.transpose() * M * w;

	Eigen::SparseMatrix<double> I(Q.rows(), Q.cols());
	I.setIdentity();

	Q = Q + 1e-8 * I;
	std::cout << "QP Solver Start: " << Q.norm() << ", " << B.norm() << ", " << Aeq.norm() << ", " << Beq.norm() << ", " << Aieq.norm() << ", " << Bieq.norm() << std::endl;
	EigenNASOQSparse qp;
	qp.setAccThresh(1e-6);
	
	delta_x.setZero();
	double perturb = 1e-9;
	bool isQPSuccess = qp.solve(Q, B, Aeq, Beq, Aieq, Bieq, lx, ux, delta_x, perturb);

	dphi = delta_x;
	double sqerrorResidual = (edgeDphi2FaceDphi * dphi - w).transpose() * M * (edgeDphi2FaceDphi * dphi - w);
	double errorResidual = std::sqrt(sqerrorResidual);
	double constraintsResidual = (C * dphi).norm();
	double initialwNorm = w.transpose() * M * w;

	std::cout << "|| w ||_ainv = " << initialwNorm << std::endl;
	std::cout << "|| dphi - w ||_ainv = " << errorResidual << std::endl;
	std::cout << "|| C * dphi || = " << constraintsResidual << std::endl;

}

void estimateWrinkleVariablesFromStrainCutbyTension(
	const std::vector<Eigen::Matrix2d>& abars,
	const Eigen::MatrixXd& curPos,
	const Eigen::MatrixXi& F,
	const std::set<int>& clampedVerts,
	double amplitudeEstimate,
	Eigen::VectorXd& amp,
	Eigen::VectorXd& phi,
	Eigen::VectorXd& dphi,
	std::set<int>& tensionFaces)
{
	tensionFaces.clear();
	MeshConnectivity mesh(F);
	MeshGeometry curGeo(curPos, mesh);
	assert(abars.size() == F.rows());
	int nfaces = F.rows();
	Eigen::MatrixXd w;
	estimateAmpOmegaFromStrain(abars, curPos, F, clampedVerts, amplitudeEstimate, amp, w, tensionFaces);
	faceDPhi2EdgeDPhi(w, tensionFaces, abars, F, dphi);

	Eigen::MatrixXd cutV;
	Eigen::MatrixXi cutF;

	Eigen::VectorXd cutAmp, cutPhi;

	roundPhiFromDphiCutbyTension(curPos, F, cutV, cutF, abars, amp, dphi, ComisoRound, phi, cutPhi, cutAmp, tensionFaces, false);// by setting to false, round phi codes will use provided pure tension faces to cut the mesh and do rounding.
	return;
}



static void integrabilityConstraintMatrix(const MeshConnectivity& mesh,
	const IntrinsicGeometry& geo,
	Eigen::SparseMatrix<double>& C)
{
	int rows = 0;
	std::vector<Eigen::Triplet<double> > Ccoeffs;
	for (int i = 0; i < mesh.nEdges(); i++)
	{
		if (mesh.edgeFace(i, 0) != -1 && mesh.edgeFace(i, 1) != -1)
		{
			int face1 = mesh.edgeFace(i, 0);
			int face2 = mesh.edgeFace(i, 1);

			int vert1 = mesh.edgeVertex(i, 0);
			int vert2 = mesh.edgeVertex(i, 1);
			
			Eigen::Vector2d barys[3] = { {0,0}, {1,0}, {0,1} };
			Eigen::Vector2d face1e(0, 0);
			for (int j = 0; j < 3; j++)
			{
				if (vert1 == mesh.faceVertex(face1, j))
					face1e -= barys[j];
				else if (vert2 == mesh.faceVertex(face1, j))
					face1e += barys[j];
			}

			face1e /= std::sqrt(face1e.transpose() * geo.abars[face1] * face1e);

			Eigen::Matrix2d T = geo.Ts.block<2, 2>(2 * i, 2);
			Eigen::RowVector2d oppdot = -face1e.transpose() * geo.abars[face1] * T;
			Eigen::RowVector2d samedot = face1e.transpose() * geo.abars[face1];
			for (int k = 0; k < 2; k++)
			{
				Ccoeffs.push_back({ rows, 2 * face1 + k, samedot[k] });
				Ccoeffs.push_back({ rows, 2 * face2 + k, oppdot[k] });
			}

			rows++;
		}
	}
	C.resize(rows, 2 * mesh.nFaces());
	C.setFromTriplets(Ccoeffs.begin(), Ccoeffs.end());
}

static void scalarDifferentialMatrix(const MeshConnectivity& mesh, const IntrinsicGeometry& geo, 
	const std::vector<int>& facesmap,
	int numsvars,
	Eigen::SparseMatrix<double>& sD)
{
	std::vector<Eigen::Triplet<double> > Dcoeffs;
	int rows = 0;
	for (int i = 0; i < mesh.nEdges(); i++)
	{
		if (mesh.edgeFace(i, 0) != -1 && mesh.edgeFace(i, 1) != -1)
		{
			int face1 = mesh.edgeFace(i, 0);
			int face2 = mesh.edgeFace(i, 1);
			if (facesmap[face1] != -1 && facesmap[face2] != -1)
			{
				Dcoeffs.push_back({ rows, facesmap[face1], 1.0 });
				Dcoeffs.push_back({ rows, facesmap[face2], -1.0 });
			}
			rows++;
		}
	}
	sD.resize(rows, numsvars);
	sD.setFromTriplets(Dcoeffs.begin(), Dcoeffs.end());
}

static void dualDifferentialMatrix(const MeshConnectivity& mesh, const IntrinsicGeometry& geo, Eigen::SparseMatrix<double>& D)
{
	std::vector<Eigen::Triplet<double> > Dcoeffs;
	int rows = 0;
	for (int i = 0; i < mesh.nEdges(); i++)
	{
		if (mesh.edgeFace(i, 0) != -1 && mesh.edgeFace(i, 1) != -1)
		{
			int face1 = mesh.edgeFace(i, 0);
			int face2 = mesh.edgeFace(i, 1);
			
			Eigen::Matrix2d T21 = geo.Ts.block<2, 2>(2 * i, 2);
			Eigen::Matrix2d T12 = geo.Ts.block<2, 2>(2 * i, 0);
			
			for (int k = 0; k < 2; k++)
			{
				Dcoeffs.push_back({ rows + k, 2 * face1 + k, 1.0 });
				for (int l = 0; l < 2; l++)
				{
					Dcoeffs.push_back({ rows + k, 2 * face2 + l, -T21(k,l) });
				}                
			}
			rows+=2;
			for (int k = 0; k < 2; k++)
			{
				Dcoeffs.push_back({ rows + k, 2 * face2 + k, 1.0 });
				for (int l = 0; l < 2; l++)
				{
					Dcoeffs.push_back({ rows + k, 2 * face1 + l, -T12(k,l) });
				}
			}
			rows += 2;
		}
	}
	D.resize(rows, 2 * mesh.nFaces());
	D.setFromTriplets(Dcoeffs.begin(), Dcoeffs.end());
}

static void dualMassMatrix(const MeshConnectivity& mesh, 
	const IntrinsicGeometry& geo, 
	Eigen::SparseMatrix<double>& M,
	const std::vector<int> &facesmap,
	double usedreg, double unusedreg)
{
	std::vector<Eigen::Triplet<double> > Mcoeffs;
	int rows = 0;
	for (int i = 0; i < mesh.nEdges(); i++)
	{
		if (mesh.edgeFace(i, 0) != -1 && mesh.edgeFace(i, 1) != -1)
		{
			int face1 = mesh.edgeFace(i, 0);
			int face2 = mesh.edgeFace(i, 1);

			double reg = unusedreg;
			if (facesmap[face1] != -1 && facesmap[face2] != -1)
				reg = usedreg;

			for (int k = 0; k < 2; k++)
				for (int l = 0; l < 2; l++)
				{
					Mcoeffs.push_back({ rows + k, rows + l, reg * geo.abars[face1](k,l) });
				}
			rows += 2;
			for (int k = 0; k < 2; k++)
				for (int l = 0; l < 2; l++)
				{
					Mcoeffs.push_back({ rows + k, rows + l, reg * geo.abars[face2](k,l) });
				}
			rows += 2;
		}
	}
	M.resize(rows, rows);
	M.setFromTriplets(Mcoeffs.begin(), Mcoeffs.end());
}

void laplacianInterpolation(const Eigen::SparseMatrix<double> L, const Eigen::VectorXd& x, const std::set<int>& clampedDOFs, Eigen::VectorXd& smoothedX)
{
	std::cout << "Laplacian interpolation." << std::endl;
	int nVars = x.size();

	Eigen::SparseMatrix<double> Proj;
	int row = 0;
	std::vector<Eigen::Triplet<double>> T;

	for (int i = 0; i < nVars; i++)
	{
		if (clampedDOFs.find(i) != clampedDOFs.end())
		{
			T.push_back(Eigen::Triplet<double>(row, i, 1.0));
			row++;
		}
	}
	Proj.resize(row, nVars);
	Proj.setFromTriplets(T.begin(), T.end());

	Eigen::VectorXd B, Beq, Bieq;
	Eigen::SparseMatrix<double> Q, Aeq, Aieq, C;
	Eigen::VectorXd lx, ux, delta_x;

	Eigen::SparseMatrix<double> I(nVars, nVars);
	I.setIdentity();

	Q = -L + 1e-8 * I;
	int nEq = row;
	int nIneq = 0;

	Aieq.resize(nIneq, nVars);
	Bieq.resize(nIneq);

	delta_x.resize(nVars);
	delta_x = x;

	Aeq = Proj;
	Beq = Proj * x;

	B.resize(nVars);
	B.setZero();


	lx.resize(nVars);
	lx.setZero();


	ux.resize(0);

	std::cout << "QP Solver Start: " << std::endl;
	EigenNASOQSparse qp;
	//if (this->m_current.fDelta > 1e-6)
	qp.setAccThresh(1e-6);
	double perturb = 1e-9;
	bool isQPSuccess = qp.solve(Q, B, Aeq, Beq, Aieq, Bieq, lx, ux, delta_x, perturb);

	smoothedX = delta_x;
	for (auto& v : clampedDOFs)
		smoothedX(v) = x(v);
}

void amplitudeSmoothing(const Eigen::VectorXd& amp, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const std::vector<int>& fxiedFaces, Eigen::VectorXd &smoothedAmp)
{
	std::set<int> freeVerts;
	for (int i = 0; i < fxiedFaces.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int vid = F(fxiedFaces[i], j);
			if (freeVerts.find(vid) == freeVerts.end())
			{
				freeVerts.insert(vid);
			}
		}
	}
	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);
	laplacianInterpolation(L, amp, freeVerts, smoothedAmp);
}


void roundPhiFromDphiCutbyTension(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	Eigen::MatrixXd& seamedV,
	Eigen::MatrixXi& seamedF,
	const std::vector<Eigen::Matrix2d>& abars,
	const Eigen::VectorXd& amp,
	const Eigen::VectorXd& dPhi, // |E| vector of phi jumps	
	const RoundingType& roundType,// how to solve the MIP, either use comiso or gurobi
	Eigen::VectorXd& phi,
	Eigen::VectorXd& seamedPhi,
	Eigen::VectorXd& seamedAmp,
	std::set<int> &problemFaces,
	bool isRecomputeProbF // The assumption is that pure tension faces having zero amplitude
	)			
{
	MeshConnectivity mesh(F);
	int nonzeroamps = 0;

	if (isRecomputeProbF || problemFaces.size() == 0)
	{
		std::cout << "recompute the problem faces. " << std::endl;
		problemFaces.clear();
		locatePotentialPureTensionFaces(abars, V, F, problemFaces);	// by default, pure tension faces have zero amplitudes


		for (int i = 0; i < F.rows(); i++)
		{
			bool isZeroFace = true;
			for (int j = 0; j < 3; j++)
			{
				if (amp(F(i, j)) > 0)
					isZeroFace = false;
			}
			if (isZeroFace)
				problemFaces.insert(i);
		}
	}

	nonzeroamps = F.rows() - problemFaces.size();

	std::vector<int> nonzero2old(nonzeroamps);
	Eigen::MatrixXi nonzeroF(nonzeroamps, 3);
	Eigen::MatrixXi zeroF(F.rows() - nonzeroamps, 3);
	std::vector<int> zero2old(F.rows() - nonzeroamps);

	int nzidx = 0;
	int zidx = 0;
	for (int i = 0; i < F.rows(); i++)
	{
		if (problemFaces.find(i) == problemFaces.end())
		{
			nonzeroF.row(nzidx) = F.row(i);
			nonzero2old[nzidx] = i;
			nzidx++;
		}
		else
		{
			zeroF.row(zidx) = F.row(i);
			zero2old[zidx] = i;
			zidx++;
		}
	}

	MeshConnectivity origNonZeroMesh(nonzeroF);

	std::vector<int> bdryedgecount(V.rows());
	for (int i = 0; i < origNonZeroMesh.nEdges(); i++)
	{
		int face1 = origNonZeroMesh.edgeFace(i, 0);
		int face2 = origNonZeroMesh.edgeFace(i, 1);
		if (face1 == -1 || face2 == -1)
		{
			bdryedgecount[origNonZeroMesh.edgeVertex(i, 0)]++;
			bdryedgecount[origNonZeroMesh.edgeVertex(i, 1)]++;
		}
	}

	std::map<int, std::vector<std::pair<int, int> > > nonmanifoldVerts;
	for (int i = 0; i < nonzeroF.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (bdryedgecount[nonzeroF(i, j)] > 2)
			{
				nonmanifoldVerts[nonzeroF(i, j)].push_back({ i,j });
			}
		}
	}

	std::vector<std::vector<std::vector<std::pair<int, int> > > > nonmanifoldClusters;
	for (auto& it : nonmanifoldVerts)
	{
		const auto& cluster = it.second;
		nonmanifoldClusters.push_back(std::vector<std::vector<std::pair<int, int> > >());
		int nclusterfaces = cluster.size();
		std::map<int, int> face2id;
		for (int i = 0; i < nclusterfaces; i++)
		{
			face2id[cluster[i].first] = i;
		}
		std::vector<bool> visited(nclusterfaces);
		for (int i = 0; i < nclusterfaces; i++)
		{
			if (visited[i])
				continue;
			std::deque<int> tovisit;
			tovisit.push_back(i);
			std::vector<std::pair<int, int> > newcluster;
			while (!tovisit.empty())
			{
				int next = tovisit.front();
				tovisit.pop_front();
				if (visited[next])
					continue;
				visited[next] = true;
				newcluster.push_back(cluster[next]);
				for (int j = 0; j < 3; j++)
				{
					int edge = origNonZeroMesh.faceEdge(cluster[next].first, j);
					int orient = origNonZeroMesh.faceEdgeOrientation(cluster[next].first, j);
					int opp = origNonZeroMesh.edgeFace(edge, 1 - orient);
					if (opp != -1)
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
			nonmanifoldClusters.back().push_back(newcluster);
		}
	}

	Eigen::MatrixXi origNonZeroF = nonzeroF;

	int freeidx = V.rows();
	for (auto& it : nonmanifoldClusters)
	{
		for (int i = 1; i < it.size(); i++)
		{
			for (auto vf : it[i])
			{
				nonzeroF(vf.first, vf.second) = freeidx;
			}
			freeidx++;
		}		
	}

	std::cout << "Found " << nonmanifoldVerts.size() << " nonmanifold vertices " << std::endl;

	MeshConnectivity punctmesh(nonzeroF);

	std::vector<std::vector<int> > cuts;
	findCuts(nonzeroF, cuts);
	std::cout << "Used " << cuts.size() << " cuts" << std::endl;

	std::map<std::pair<int, int>, int> cutpairs;

	// convert cuts to edge indices
	for (int i = 0; i < cuts.size(); i++)
	{
		int len = cuts[i].size();
		for (int j = 0; j < len - 1; j++)
		{
			cutpairs[std::pair<int, int>(cuts[i][j], cuts[i][j + 1])] = i;
			//std::cout << cuts[i][j] << ", ";
		}
	}

	Eigen::MatrixXi newF;

	cutMesh(nonzeroF, cuts, newF);
	reindex(newF);

	int newverts = 0;
	for (int i = 0; i < newF.rows(); i++)
		for (int j = 0; j < 3; j++)
			newverts = std::max(newverts, newF(i, j) + 1);

	reindex(zeroF);
	int zeroverts = 0;
	for (int i = 0; i < zeroF.rows(); i++)
		for (int j = 0; j < 3; j++)
			zeroverts = std::max(zeroverts, zeroF(i, j) + 1);

	int intdofs = 0;	

	std::vector<Eigen::Triplet<double> > Ccoeffs;

	// integer constraints on the two sides of the cut
	int row = 0;
	for (int i = 0; i < punctmesh.nEdges(); i++)
	{
		int f1 = punctmesh.edgeFace(i, 0);
		int f2 = punctmesh.edgeFace(i, 1);
		if (f1 == -1 || f2 == -1)
			continue;
		int v1 = punctmesh.edgeVertex(i, 0);
		int v2 = punctmesh.edgeVertex(i, 1);
		double sign = 1.0;
		auto it = cutpairs.find(std::pair<int, int>(v1, v2));
		if (it == cutpairs.end())
		{
			it = cutpairs.find(std::pair<int, int>(v2, v1));
			sign = -1.0;
		}
		if (it == cutpairs.end())
			continue;

		int newv1 = -1;
		int newv2 = -1;
		int neww1 = -1;
		int neww2 = -1;
		for (int j = 0; j < 3; j++)
		{
			if (nonzeroF(f1, j) == v1)
				newv1 = newF(f1, j);
			if (nonzeroF(f2, j) == v1)
				newv2 = newF(f2, j);
			if (nonzeroF(f1, j) == v2)
				neww1 = newF(f1, j);
			if (nonzeroF(f2, j) == v2)
				neww2 = newF(f2, j);
		}

		Ccoeffs.push_back(Eigen::Triplet<double>(row, newv1, sign));
		Ccoeffs.push_back(Eigen::Triplet<double>(row, newv2, -sign));
		Ccoeffs.push_back(Eigen::Triplet<double>(row, newverts + intdofs, 1.0));
		Ccoeffs.push_back(Eigen::Triplet<double>(row + 1, neww1, sign));
		Ccoeffs.push_back(Eigen::Triplet<double>(row + 1, neww2, -sign));
		Ccoeffs.push_back(Eigen::Triplet<double>(row + 1, newverts + intdofs, 1.0));
		row += 2;
		intdofs++;
	}

	// add integer jumps for non-manifold vertices that connect the same connected component to itself

	std::vector<int> parent(newverts);
	for (int i = 0; i < newverts; i++)
		parent[i] = i;

	for (int i = 0; i < newF.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int idx1 = newF(i, j);
			int idx2 = newF(i, (j + 1) % 3);
			while (parent[idx1] != idx1)
				idx1 = parent[idx1];
			while (parent[idx2] != idx2)
				idx2 = parent[idx2];
			if (idx1 != idx2)
			{
				// union
				parent[idx2] = idx1;
			}
		}
	}

	int redundantManiJumps = 0;
	int essentialManiJumps = 0;

	for (auto &it : nonmanifoldClusters)
	{
		int nclusters = it.size();		
		assert(nclusters > 1);
		int basef = it[0][0].first;
		int basev = it[0][0].second;

		for (int j = 1; j < nclusters; j++)
		{
			int newv1 = newF(basef, basev);
			int newv2 = newF(it[j][0].first, it[j][0].second);

			int idx1 = newv1;
			int idx2 = newv2;
			while (parent[idx1] != idx1)
				idx1 = parent[idx1];
			while (parent[idx2] != idx2)
				idx2 = parent[idx2];
			if (idx1 == idx2)
			{
				essentialManiJumps++;
				Ccoeffs.push_back({ row, newv1, 1.0 });
				Ccoeffs.push_back({ row, newv2, -1.0 });
				Ccoeffs.push_back({ row, newverts + intdofs, 1.0 });
				row++;
				intdofs++;
			}
			else
			{
				// integer jump not needed; force 0
				redundantManiJumps++;
				parent[idx2] = idx1;
				Ccoeffs.push_back({ row, newv1, 1.0 });
				Ccoeffs.push_back({ row, newv2, -1.0 });
				row++;
			}
		}
	}
	Eigen::SparseMatrix<double> C(row, newverts + intdofs + 1);
	C.setFromTriplets(Ccoeffs.begin(), Ccoeffs.end());

	std::cout << "For non-manifold vertices, added " << essentialManiJumps << " integer jumps and " << redundantManiJumps << " zero-jumps" << std::endl;

	std::vector<Eigen::Triplet<double> > Acoeffs;
	int newfaces = newF.rows();
	//assert(newfaces == punctF.rows());

	for (int i = 0; i < newfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Acoeffs.push_back(Eigen::Triplet<double>(3 * i + j, newF(i, j), 1.0));
		}
	}
	Eigen::SparseMatrix<double> A(3 * newfaces, newverts + intdofs);
	A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

	std::vector<Eigen::Triplet<double> > L2Minvcoeffs;
	for (int i = 0; i < newfaces; i++)
	{
		double area = 0.5 * std::sqrt(abars[nonzero2old[i]].determinant());
		for (int j = 0; j < 3; j++)
		{
			//double avamp = 0.5 * (amp[F(nonzero2old[i], (j + 1) % 3)] + amp[F(nonzero2old[i], (j + 2) % 3)]);
			//L2Minvcoeffs.push_back(Eigen::Triplet<double>(3 * i + j, 3 * i + j, area * avamp));
			L2Minvcoeffs.push_back(Eigen::Triplet<double>(3 * i + j, 3 * i + j, area));
		}
	}
	Eigen::SparseMatrix<double> L2Minv(3 * newfaces, 3 * newfaces);
	L2Minv.setFromTriplets(L2Minvcoeffs.begin(), L2Minvcoeffs.end());



	std::vector<Eigen::Triplet<double> > Dcoeffs;
	for (int i = 0; i < newfaces; i++)
	{
		Dcoeffs.push_back(Eigen::Triplet<double>(3 * i, 3 * i + 2, 1.0));
		Dcoeffs.push_back(Eigen::Triplet<double>(3 * i, 3 * i + 1, -1.0));
		Dcoeffs.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i, 1.0));
		Dcoeffs.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 2, -1.0));
		Dcoeffs.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 1, 1.0));
		Dcoeffs.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i, -1.0));
	}
	Eigen::SparseMatrix<double> D(3 * newfaces, 3 * newfaces);
	D.setFromTriplets(Dcoeffs.begin(), Dcoeffs.end());

	Eigen::SparseMatrix<double> Mat = A.transpose() * D.transpose() * L2Minv * D * A;

	Eigen::VectorXd punctfield(3 * newfaces);
	for (int i = 0; i < newfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int edgeidx = mesh.faceEdge(nonzero2old[i], j);
			double sign = (mesh.faceEdgeOrientation(nonzero2old[i], j) == 0 ? 1.0 : -1.0);
			punctfield[3 * i + j] = dPhi[edgeidx] * sign / 2.0 / PI;
		}
	}
	Eigen::VectorXd rhs = A.transpose() * D.transpose() * L2Minv * punctfield;

	Eigen::VectorXd result;
	Eigen::VectorXi toRound(intdofs);
	for (int i = 0; i < intdofs; i++)
	{
		toRound[i] = newverts + i;
	}
	std::cout << intdofs << " integer variables" << std::endl;
	std::cout << "Comiso round" << std::endl;
	ComisoWrapper(C, Mat, result, rhs, toRound, 1e-6);
		
	std::cout << "Solver residual: " << (Mat * result - rhs).norm() << std::endl;
	std::cout << "Reconstruction residual: " << (D * A * result - punctfield).transpose() * L2Minv * (D * A * result - punctfield) << std::endl;
	std::cout << "Integer jumps: ";
	for (int i = 0; i < intdofs; i++)
		std::cout << result[newverts + i] << " ";
	std::cout << std::endl;

	// map phi back to the original mesh
	int oldverts = V.rows();
	phi.resize(oldverts);
	phi.setZero();
	for (int i = 0; i < newfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			phi[origNonZeroF(i, j)] = 2.0 * PI * result[newF(i, j)];
		}
	}

	seamedV.resize(newverts + zeroverts, 3);

	// construct new position
	for (int i = 0; i < newfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			seamedV.row(newF(i, j)) = V.row(origNonZeroF(i, j));
		}
	}

	for (int i = 0; i < zeroF.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			seamedV.row(newverts + zeroF(i, j)) = V.row(F(zero2old[i], j));
		}
	}

	seamedF.resize(F.rows(), 3);
	for (int i = 0; i < newfaces; i++)
	{
		seamedF.row(nonzero2old[i]) = newF.row(i);
	}
	for (int i = 0; i < zeroF.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			seamedF(zero2old[i], j) = newverts + zeroF(i, j);
		}
	}

	// construct new phi and amp with the seam
	seamedPhi.resize(newverts + zeroverts);
	seamedAmp.resize(newverts + zeroverts);
	seamedPhi.setZero();
	seamedAmp.setZero();
	std::set<int> freeDOFs;
	for (int i = 0; i < newfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			seamedPhi(newF(i, j)) = 2.0 * PI * result(newF(i, j));
			seamedAmp(newF(i, j)) = amp(origNonZeroF(i, j));
		}
	}

	for (int i = 0; i < zeroF.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			//seamedAmp(newverts + zeroF(i, j)) = amp(F(zero2old[i], j));
			seamedAmp(newverts + zeroF(i, j)) = 0;
		}
	}
	//igl::writeOBJ("../cut.obj", seamedV, seamedF);
}

void roundPhiFromOmegaCutbyTension(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	Eigen::MatrixXd& seamedV,
	Eigen::MatrixXi& seamedF,
	const std::vector<Eigen::Matrix2d>& abars,
	const Eigen::VectorXd& amp,
	const Eigen::MatrixXd& w, // |F| x 2 one-form in the barycentric basis, i.e. w.row(i){1,0} is the integral of w along edge (0,1) of the ith face
	Eigen::VectorXd& phi,
	Eigen::VectorXd& seamedPhi,
	Eigen::VectorXd& seamedAmp,
	std::set<int>& problemFaces,
	bool isRecomputeProbF)
{
	MeshConnectivity mesh(F);
	int nonzeroamps = 0;

	if (isRecomputeProbF || problemFaces.size() == 0)
	{
		problemFaces.clear();
		locatePotentialPureTensionFaces(abars, V, F, problemFaces);	// by default, pure tension faces have zero amplitudes


		for (int i = 0; i < F.rows(); i++)
		{
			bool isZeroFace = true;
			for (int j = 0; j < 3; j++)
			{
				if (amp(F(i, j)) > 0)
					isZeroFace = false;
			}
			if (isZeroFace)
			{
				problemFaces.insert(i);
			}
		}
	}
	

	nonzeroamps = F.rows() - problemFaces.size();

	std::vector<int> nonzero2old(nonzeroamps);
	Eigen::MatrixXi nonzeroF(nonzeroamps, 3);
	Eigen::MatrixXi zeroF(F.rows() - nonzeroamps, 3);
	std::vector<int> zero2old(F.rows() - nonzeroamps);

	int nzidx = 0;
	int zidx = 0;
	for (int i = 0; i < F.rows(); i++)
	{
		if (problemFaces.find(i) == problemFaces.end())
		{
			nonzeroF.row(nzidx) = F.row(i);
			nonzero2old[nzidx] = i;
			nzidx++;
		}
		else
		{
			zeroF.row(zidx) = F.row(i);
			zero2old[zidx] = i;
			zidx++;
		}
	}

	MeshConnectivity origNonZeroMesh(nonzeroF);

	std::vector<int> bdryedgecount(V.rows());
	for (int i = 0; i < origNonZeroMesh.nEdges(); i++)
	{
		int face1 = origNonZeroMesh.edgeFace(i, 0);
		int face2 = origNonZeroMesh.edgeFace(i, 1);
		if (face1 == -1 || face2 == -1)
		{
			bdryedgecount[origNonZeroMesh.edgeVertex(i, 0)]++;
			bdryedgecount[origNonZeroMesh.edgeVertex(i, 1)]++;
		}
	}

	std::map<int, std::vector<std::pair<int, int> > > nonmanifoldVerts;
	for (int i = 0; i < nonzeroF.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (bdryedgecount[nonzeroF(i, j)] > 2)
			{
				nonmanifoldVerts[nonzeroF(i, j)].push_back({ i,j });
			}
		}
	}

	std::vector<std::vector<std::vector<std::pair<int, int> > > > nonmanifoldClusters;
	for (auto& it : nonmanifoldVerts)
	{
		const auto& cluster = it.second;
		nonmanifoldClusters.push_back(std::vector<std::vector<std::pair<int, int> > >());
		int nclusterfaces = cluster.size();
		std::map<int, int> face2id;
		for (int i = 0; i < nclusterfaces; i++)
		{
			face2id[cluster[i].first] = i;
		}
		std::vector<bool> visited(nclusterfaces);
		for (int i = 0; i < nclusterfaces; i++)
		{
			if (visited[i])
				continue;
			std::deque<int> tovisit;
			tovisit.push_back(i);
			std::vector<std::pair<int, int> > newcluster;
			while (!tovisit.empty())
			{
				int next = tovisit.front();
				tovisit.pop_front();
				if (visited[next])
					continue;
				visited[next] = true;
				newcluster.push_back(cluster[next]);
				for (int j = 0; j < 3; j++)
				{
					int edge = origNonZeroMesh.faceEdge(cluster[next].first, j);
					int orient = origNonZeroMesh.faceEdgeOrientation(cluster[next].first, j);
					int opp = origNonZeroMesh.edgeFace(edge, 1 - orient);
					if (opp != -1)
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
			nonmanifoldClusters.back().push_back(newcluster);
		}
	}

	Eigen::MatrixXi origNonZeroF = nonzeroF;

	int freeidx = V.rows();
	for (auto& it : nonmanifoldClusters)
	{
		for (int i = 1; i < it.size(); i++)
		{
			for (auto vf : it[i])
			{
				nonzeroF(vf.first, vf.second) = freeidx;
			}
			freeidx++;
		}
	}

	std::cout << "Found " << nonmanifoldVerts.size() << " nonmanifold vertices " << std::endl;

	MeshConnectivity punctmesh(nonzeroF);

	std::vector<std::vector<int> > cuts;
	findCuts(nonzeroF, cuts);
	std::cout << "Used " << cuts.size() << " cuts" << std::endl;

	std::map<std::pair<int, int>, int> cutpairs;

	// convert cuts to edge indices
	for (int i = 0; i < cuts.size(); i++)
	{
		int len = cuts[i].size();
		for (int j = 0; j < len - 1; j++)
		{
			cutpairs[std::pair<int, int>(cuts[i][j], cuts[i][j + 1])] = i;
			//std::cout << cuts[i][j] << ", ";
		}
	}

	Eigen::MatrixXi newF;

	cutMesh(nonzeroF, cuts, newF);
	reindex(newF);

	int newverts = 0;
	for (int i = 0; i < newF.rows(); i++)
		for (int j = 0; j < 3; j++)
			newverts = std::max(newverts, newF(i, j) + 1);

	reindex(zeroF);
	int zeroverts = 0;
	for (int i = 0; i < zeroF.rows(); i++)
		for (int j = 0; j < 3; j++)
			zeroverts = std::max(zeroverts, zeroF(i, j) + 1);

	int intdofs = 0;

	std::vector<Eigen::Triplet<double> > Ccoeffs;

	// integer constraints on the two sides of the cut
	int row = 0;
	for (int i = 0; i < punctmesh.nEdges(); i++)
	{
		int f1 = punctmesh.edgeFace(i, 0);
		int f2 = punctmesh.edgeFace(i, 1);
		if (f1 == -1 || f2 == -1)
			continue;
		int v1 = punctmesh.edgeVertex(i, 0);
		int v2 = punctmesh.edgeVertex(i, 1);
		double sign = 1.0;
		auto it = cutpairs.find(std::pair<int, int>(v1, v2));
		if (it == cutpairs.end())
		{
			it = cutpairs.find(std::pair<int, int>(v2, v1));
			sign = -1.0;
		}
		if (it == cutpairs.end())
			continue;

		int newv1 = -1;
		int newv2 = -1;
		int neww1 = -1;
		int neww2 = -1;
		for (int j = 0; j < 3; j++)
		{
			if (nonzeroF(f1, j) == v1)
				newv1 = newF(f1, j);
			if (nonzeroF(f2, j) == v1)
				newv2 = newF(f2, j);
			if (nonzeroF(f1, j) == v2)
				neww1 = newF(f1, j);
			if (nonzeroF(f2, j) == v2)
				neww2 = newF(f2, j);
		}

		Ccoeffs.push_back(Eigen::Triplet<double>(row, newv1, sign));
		Ccoeffs.push_back(Eigen::Triplet<double>(row, newv2, -sign));
		Ccoeffs.push_back(Eigen::Triplet<double>(row, newverts + intdofs, 1.0));
		Ccoeffs.push_back(Eigen::Triplet<double>(row + 1, neww1, sign));
		Ccoeffs.push_back(Eigen::Triplet<double>(row + 1, neww2, -sign));
		Ccoeffs.push_back(Eigen::Triplet<double>(row + 1, newverts + intdofs, 1.0));
		row += 2;
		intdofs++;
	}
	for (auto& it : nonmanifoldClusters)
	{
		int nclusters = it.size();
		assert(nclusters > 1);
		int basef = it[0][0].first;
		int basev = it[0][0].second;

		for (int j = 1; j < nclusters; j++)
		{
			int newv1 = newF(basef, basev);
			int newv2 = newF(it[j][0].first, it[j][0].second);
			Ccoeffs.push_back({ row, newv1, 1.0 });
			Ccoeffs.push_back({ row, newv2, -1.0 });
			Ccoeffs.push_back({ row, newverts + intdofs, 1.0 });
			row++;
			intdofs++;
		}
	}
	Eigen::SparseMatrix<double> C(row, newverts + intdofs + 1);
	C.setFromTriplets(Ccoeffs.begin(), Ccoeffs.end());

	std::vector<Eigen::Triplet<double> > Acoeffs;
	int newfaces = newF.rows();
	//assert(newfaces == punctF.rows());

	for (int i = 0; i < newfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Acoeffs.push_back(Eigen::Triplet<double>(3 * i + j, newF(i, j), 1.0));
		}
	}
	Eigen::SparseMatrix<double> A(3 * newfaces, newverts + intdofs);
	A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

	std::vector<Eigen::Triplet<double> > L2Minvcoeffs;
	for (int i = 0; i < newfaces; i++)
	{
		double area = 0.5 * std::sqrt(abars[nonzero2old[i]].determinant());
		for (int j = 0; j < 2; j++)
		{
			//double avamp = 0.5 * (amp[F(nonzero2old[i], (j + 1) % 3)] + amp[F(nonzero2old[i], (j + 2) % 3)]);
			//L2Minvcoeffs.push_back(Eigen::Triplet<double>(3 * i + j, 3 * i + j, area * avamp));
			L2Minvcoeffs.push_back(Eigen::Triplet<double>(2 * i + j, 2 * i + j, area));
		}
	}
	Eigen::SparseMatrix<double> L2Minv(2 * newfaces, 2 * newfaces);
	L2Minv.setFromTriplets(L2Minvcoeffs.begin(), L2Minvcoeffs.end());



	std::vector<Eigen::Triplet<double> > Dcoeffs;
	for (int i = 0; i < newfaces; i++)
	{
		Dcoeffs.push_back(Eigen::Triplet<double>(2 * i, 3 * i + 1, 1.0));
		Dcoeffs.push_back(Eigen::Triplet<double>(2 * i, 3 * i, -1.0));
		Dcoeffs.push_back(Eigen::Triplet<double>(2 * i + 1, 3 * i + 2, 1.0));
		Dcoeffs.push_back(Eigen::Triplet<double>(2 * i + 1, 3 * i, -1.0));
	}

	Eigen::SparseMatrix<double> D(2 * newfaces, 3 * newfaces);
	D.setFromTriplets(Dcoeffs.begin(), Dcoeffs.end());

	Eigen::SparseMatrix<double> Mat = A.transpose() * D.transpose() * L2Minv * D * A;
	Eigen::VectorXd punctfield(2 * newfaces);
	for (int i = 0; i < newfaces; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			punctfield[2 * i + j] = w(nonzero2old[i], j) / 2.0 / PI;
		}
	}

	Eigen::VectorXd rhs = A.transpose() * D.transpose() * L2Minv * punctfield;

	Eigen::VectorXd result;
	Eigen::VectorXi toRound(intdofs);
	for (int i = 0; i < intdofs; i++)
	{
		toRound[i] = newverts + i;
	}

	std::cout << intdofs << " integer variables" << std::endl;
	ComisoWrapper(C, Mat, result, rhs, toRound, 1e-6);
	std::cout << "Solver residual: " << (Mat * result - rhs).norm() << std::endl;
	std::cout << "Reconstruction residual: " << (D * A * result - punctfield).transpose() * L2Minv * (D * A * result - punctfield) << std::endl;
	std::cout << "Integer jumps: ";
	for (int i = 0; i < intdofs; i++)
		std::cout << result[newverts + i] << " ";
	std::cout << std::endl;

	// map phi back to the original mesh
	int oldverts = V.rows();
	phi.resize(oldverts);
	phi.setZero();
	for (int i = 0; i < newfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			phi[origNonZeroF(i, j)] = 2.0 * PI * result[newF(i, j)];
		}
	}

	seamedV.resize(newverts + zeroverts, 3);

	// construct new position
	for (int i = 0; i < newfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			seamedV.row(newF(i, j)) = V.row(origNonZeroF(i, j));
		}
	}

	for (int i = 0; i < zeroF.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			seamedV.row(newverts + zeroF(i, j)) = V.row(F(zero2old[i], j));
		}
	}

	seamedF.resize(F.rows(), 3);
	for (int i = 0; i < newfaces; i++)
	{
		seamedF.row(nonzero2old[i]) = newF.row(i);
	}
	for (int i = 0; i < zeroF.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			seamedF(zero2old[i], j) = newverts + zeroF(i, j);
		}
	}

	// construct new phi and amp with the seam
	seamedPhi.resize(newverts + zeroverts);
	seamedAmp.resize(newverts + zeroverts);
	seamedPhi.setZero();
	seamedAmp.setZero();
	std::set<int> freeDOFs;
	for (int i = 0; i < newfaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			seamedPhi(newF(i, j)) = 2.0 * PI * result(newF(i, j));
			seamedAmp(newF(i, j)) = amp(origNonZeroF(i, j));
		}
	}

	for (int i = 0; i < zeroF.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			//seamedAmp(newverts + zeroF(i, j)) = amp(F(zero2old[i], j));
			seamedAmp(newverts + zeroF(i, j)) = 0;
		}
	}
}