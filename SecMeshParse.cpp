#include "SecMeshParse.h"

std::map<std::pair<int, int>, int> he2Edge(const Eigen::MatrixXi& faces)
{
	std::map< std::pair<int, int>, int > heToEdge;
	std::vector< std::vector<int> > edgeToVert;
	for (int face = 0; face < faces.rows(); ++face)
	{
		for (int i = 0; i < 3; ++i)
		{
			int vi = faces(face, i);
			int vj = faces(face, (i + 1) % 3);
			assert(vi != vj);

			std::pair<int, int> he = std::make_pair(vi, vj);
			if (he.first > he.second) std::swap(he.first, he.second);
			if (heToEdge.find(he) != heToEdge.end()) continue;

			heToEdge[he] = edgeToVert.size();
			edgeToVert.push_back(std::vector<int>(2));
			edgeToVert.back()[0] = he.first;
			edgeToVert.back()[1] = he.second;
		}
	}
	return heToEdge;
}

std::map<std::pair<int, int>, int> he2Edge(const std::vector< std::vector<int>>& edgeToVert)
{
	std::map< std::pair<int, int>, int > heToEdge;
	for (int i = 0; i < edgeToVert.size(); i++)
	{
		std::pair<int, int> he = std::make_pair(edgeToVert[i][0], edgeToVert[i][1]);
		heToEdge[he] = i;
	}
	return heToEdge;
}

Eigen::VectorXd swapEdgeVec(const std::vector< std::vector<int>>& edgeToVert, const Eigen::VectorXd& edgeVec, int flag)
{
	Eigen::VectorXd  edgeVecSwap = edgeVec;
	std::map< std::pair<int, int>, int > heToEdge = he2Edge(edgeToVert);

	int idx = 0;
	for (auto it : heToEdge)
	{
		if (flag == 0)   // ours to secstencils
			edgeVecSwap(it.second) = edgeVec(idx);
		else
			edgeVecSwap(idx) = edgeVec(it.second);
		idx++;
	}
	return edgeVecSwap;
}

Eigen::VectorXd swapEdgeVec(const Eigen::MatrixXi& faces, const Eigen::VectorXd& edgeVec, int flag)
{
	Eigen::VectorXd  edgeVecSwap = edgeVec;
	std::map< std::pair<int, int>, int > heToEdge = he2Edge(faces);

	int idx = 0;
	for (auto it : heToEdge)
	{
		if (flag == 0)   // ours to secstencils
			edgeVecSwap(it.second) = edgeVec(idx);
		else
			edgeVecSwap(idx) = edgeVec(it.second);
		idx++;
	}
	return edgeVecSwap;
}

std::vector<std::vector<int>> swapEdgeIndices(const Eigen::MatrixXi& faces, const std::vector<std::vector<int>>& edgeIndices, int flag)
{
	std::vector<std::vector<int>> edgeIndicesSwap = edgeIndices;
	std::map< std::pair<int, int>, int > heToEdge = he2Edge(faces);

	int idx = 0;
	for (auto it : heToEdge)
	{
		if (flag == 0)   // ours to secstencils
		{
			edgeIndicesSwap[it.second] = edgeIndices[idx];
		}
		else
		{
			edgeIndicesSwap[idx] = edgeIndices[it.second];
		}
		idx++;
	}

	return edgeIndicesSwap;
}

Eigen::MatrixXd edgeVec2FaceVec(const Mesh& mesh, Eigen::VectorXd& edgeVec)
{
	int nfaces = mesh.GetFaceCount();
	int nedges = mesh.GetEdgeCount();
	Eigen::MatrixXd fVec(nfaces, 3);
	fVec.setZero();

	for (int f = 0; f < nfaces; f++)
	{
		std::vector<int> faceEdges = mesh.GetFaceEdges(f);
		std::vector<int> faceVerts = mesh.GetFaceVerts(f);
		for (int j = 0; j < 3; j++)
		{
			int vid = faceVerts[j];
			int eid0 = faceEdges[j];
			int eid1 = faceEdges[(j + 2) % 3];

			Eigen::Vector3d e0 = mesh.GetVertPos(faceVerts[(j + 1) % 3]) - mesh.GetVertPos(vid);
			Eigen::Vector3d e1 = mesh.GetVertPos(faceVerts[(j + 2) % 3]) - mesh.GetVertPos(vid);

			int flag0 = 1, flag1 = 1;
			Eigen::Vector2d rhs;

			if (mesh.GetEdgeVerts(eid0)[0] == vid)
			{
				flag0 = 1;
			}
			else
			{
				flag0 = -1;
			}


			if (mesh.GetEdgeVerts(eid1)[0] == vid)
			{
				flag1 = 1;
			}
			else
			{
				flag1 = -1;
			}
			rhs(0) = flag0 * edgeVec(eid0);
			rhs(1) = flag1 * edgeVec(eid1);

			Eigen::Matrix2d I;
			I << e0.dot(e0), e0.dot(e1), e1.dot(e0), e1.dot(e1);
			Eigen::Vector2d sol = I.inverse() * rhs;

			fVec.row(f) += (sol(0) * e0 + sol(1) * e1) / 3;
		}
	}
	return fVec;
}

Mesh convert2SecMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
	std::vector<Eigen::Vector3d> pos;
	std::vector<std::vector<int>> faces;

	pos.resize(V.rows());
	for (int i = 0; i < V.rows(); i++)
	{
		pos[i] = V.row(i);
	}

	faces.resize(F.rows());
	for (int i = 0; i < F.rows(); i++)
	{
		faces[i] = { F(i, 0), F(i, 1), F(i, 2) };
	}

	Mesh mesh;

	mesh.Populate(pos, faces);
	return mesh;
}

void parseSecMesh(const Mesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
	mesh.GetPos(V);
	int nfaces = mesh.GetFaceCount();
	F.resize(nfaces, 3);

	for (int i = 0; i < nfaces; i++)
	{
		F.row(i) << mesh.GetFaceVerts(i)[0], mesh.GetFaceVerts(i)[1], mesh.GetFaceVerts(i)[2];
	}
}