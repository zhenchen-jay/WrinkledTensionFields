#include "MeshConnectivity.h"
#include <igl/boundary_facets.h>
#include <vector>
#include <stack>
#include <map>

MeshConnectivity::MeshConnectivity()
{
    F.resize(0, 3);
    FE.resize(0, 3);
    FEorient.resize(0, 3);
    EV.resize(0, 2);
    EF.resize(0, 2);
    EOpp.resize(0, 2);
}

MeshConnectivity::MeshConnectivity(const Eigen::MatrixXi &F) : F(F)
{
    std::map<std::pair<int, int>, Eigen::Vector2i > edgeFaces;
    int nfaces = F.rows();
    for(int i=0; i<nfaces; i++)
    {
        for(int j=0; j<3; j++)
        {
            int v0 = F(i, (j+1)%3);
            int v1 = F(i, (j+2)%3);
            int idx=0;
            if(v0 > v1) 
            {
                std::swap(v0,v1);
                idx = 1;
            }
            
            std::pair<int, int> p(v0,v1);
            auto it = edgeFaces.find(p);
            if(it == edgeFaces.end())
            {
                edgeFaces[p][idx] = i;
                edgeFaces[p][1-idx] = -1;
            }
            else
            {
                edgeFaces[p][idx] = i;
            }
        }
    }
    
    int nedges = edgeFaces.size();
    FE.resize(nfaces, 3);
    FEorient.resize(nfaces, 3);
    EV.resize(nedges, 2);
    EF.resize(nedges, 2);
    EOpp.resize(nedges, 2);
    std::map<std::pair<int, int>, int> edgeIndices;
    
    int idx=0;
    for(auto it : edgeFaces)
    {
        edgeIndices[it.first] = idx;
        EV(idx, 0) = it.first.first;
        EV(idx, 1) = it.first.second;
        EF(idx, 0) = it.second[0];
        EF(idx, 1) = it.second[1];
        idx++;
    }
    
    Eigen::MatrixXi bnd;
    igl::boundary_facets(F, bnd);
    for(int i=0; i<nfaces; i++)
    {
        for(int j=0; j<3; j++)
        {
            int v0 = F(i, (j+1)%3);
            int v1 = F(i, (j+2)%3);
            if(v0 > v1) std::swap(v0,v1);
            FE(i,j) = edgeIndices[std::pair<int,int>(v0,v1)];
            for (int k = 0; k < bnd.rows(); k++)
            {
                int ev0 = bnd(k,0);
                int ev1 = bnd(k,1);
                if (ev0 > ev1) std::swap(ev0,ev1);
                if (ev0 == v0 && ev1 == v1)
                {
                    bnd_edges.push_back(Eigen::Vector3i(bnd(k,0), bnd(k,1), F(i,(j+3)%3)));
                }
            }
        }
    }

    
    for (int i = 0; i < nedges; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            EOpp(i, j) = oppositeVertex(i, j);
        }
    }

    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int edge = faceEdge(i, j);
            if (edgeFace(edge, 0) == i)
                FEorient(i, j) = 0;
            else
                FEorient(i, j) = 1;
        }
    }
}

int MeshConnectivity::oppositeVertexIndex(int edge, int faceidx) const
{
    int face = edgeFace(edge, faceidx);
    if(face == -1)
        return -1;
        
    for(int j=0; j<3; j++)
    {
        if(F(face, j) != edgeVertex(edge, 0) && F(face, j) != edgeVertex(edge, 1))
            return j;
    }    
    return -1;
}

int MeshConnectivity::oppositeVertex(int edge, int faceidx) const
{
    int face = edgeFace(edge, faceidx);
    int idx = oppositeVertexIndex(edge, faceidx);
    if(idx == -1)
        return -1;
    return F(face, idx);
}

int MeshConnectivity::vertexOppositeFaceEdge(int face, int vertidx) const
{
    int edge = faceEdge(face, vertidx);
    int edgeorient = faceEdgeOrientation(face, vertidx);
    return edgeOppositeVertex(edge, 1 - edgeorient);
}

/*
bool MeshConnectivity::constructMST(int nverts)
{
    int nedges = nEdges();
    edgeParentID.resize(nedges);
    std::vector<bool> visitedV(nverts,false);
    std::vector<bool> visitedEdge(nedges,false);
    std::queue<int> edgeList;
    std::vector<int> MST;

    visitedEdge[0] = true;
    visitedV[edgeVertex(0,0)] = true;
    visitedV[edgeVertex(0,1)] = true;
    edgeList.push(0);
    MST.push_back(0);
    edgeParentID[0].push_back(0);
    
    do
    {
        int curEID = edgeList.top();
        edgeList.pop();
        visitedEdge[curEID] = true;
        for (int i = 0; i < 2; i++)
        {
            int faceID = edgeFace(curEID,i);
            if (faceID >= 0)
            {
                for (int j = 0; j < 3; j++)
                {
                    int eid = faceEdge(faceID,j);
                    if (!visitedEdge[eid])
                        edgeList.push(eid);
                }
            }
        }
        //if (!visitedEdge[curEID])
        //{
            if (visitedV[edgeVertex(curEID,0)] && (!visitedV[edgeVertex(curEID,1)]))
            {
                MST.push_back(curEID);
                visitedV[edgeVertex(curEID,1)] = true;
            //    edgeParentID[curEID].push_back(curEID); 
            }
            else if ((!visitedV[edgeVertex(curEID,0)]) && visitedV[edgeVertex(curEID,1)])
            {
                MST.push_back(curEID);
                visitedV[edgeVertex(curEID,0)] = true;
            //    edgeParentID[curEID].push_back(curEID); 
            }
            //else if (visitedV[edgeVertex(curEID,0)] && visitedV[edgeVertex(curEID,1)])
            //{
            //    int faceID = -1;
            //    int vopp = -1;
            //    //find parent ID
            //    for (int j = 0; j < 2; j++)
            //    {
            //        int vopp_temp = edgeOppositeVertex(curEID,j);
            //        if (vopp_temp >= 0)
            //        {
            //            if (visitedV[vopp_temp])
            //            {
            //                faceID = edgeFace(curEID,j);
            //                vopp = vopp_temp;
            //            }
            //        }
            //    }
            //    if (faceID < 0 || vopp < 0)
            //    {
            //        std::cout << "minimum spanning tree construction failed!" << std::endl;
            //        return false;
            //    }
            //    for (int j = 0; j < 3 ; j++)
            //    {
            //        int parentEdge = faceEdge(faceID, j);
            //        //std::cout << "parent edge vertex : " << edgeVertex(parentEdge,0) << " " << edgeVertex(parentEdge,1) << std::endl;
            //        if (vopp == edgeVertex(parentEdge,0) || vopp == edgeVertex(parentEdge,1) )
            //        {
            //            edgeParentID[curEID].push_back(parentEdge); 
            //        }
            //    }
            //}
            else if ((!visitedV[edgeVertex(curEID,0)]) && (!visitedV[edgeVertex(curEID,1)]))
            {
                std::cout << "Something is wrong when constructing the minimum spanning tree!!!!" << std::endl;
                return false;    
            } 
        //}
        //visitedEdge[curEID] = true;
    }while(!edgeList.empty());

    if (MST.size() != nverts -1)
    {
        std::cout << "MST construction is wrong!" << std::endl;
        return false;
    }

    //Construct the reduce basis
    Eigen::MatrixXi reducedB(nedges, MST.size());
    reducedB.setZero();
    std::vector<bool> assignedEdges(nedges,false);
    std::sort(MST.begin(), MST.end());

    std::vector<int> remainingEdges;
    
    for (int i = 0; i < nedges; i++)
    {
        int cnt = 0;
        if (i == MST[cnt])
        {
            reducedB(i,cnt) = 1;    
            assignedEdges[i] = true;
            cnt ++;
        }
        else
            remainingEdges.push_back(i); 
    }

    while (remainingEdges.size() > 0)
    {
        for (int i = 0; i < remainingEdges.size(); i++)
        {
            int curEid = remainingEdges[i];
            for (int k = 0; k < 2; k++)
            {
                int fid = edgeFace(curEid,k);
                int e[2];
                int cnt = 0;
                for (int j = 0; j < 3; j++)
                {
                    int edge = faceEdge(fid,j);
                    if (edge!=curEid)
                    {
                        e[cnt] = edge;
                        cnt ++;
                    }
                }
                if (assignedEdge[e[0]] == true && assignedEdge[e[1]] == true)
                {
                    reducedB.row(curEid) = -(reducedB.row(e[0]) + reducedB.row(e[1]));
                    assignedEdge[curEid] = true;
                    //remove from list
                    std::swap(remainingEdge.begin() + i,remainingEdges.back());
                    remainingEdges.pop_back();
                    i--;
                    break;
                }
            }
        }
    }
    



    //std::cout << "number of verts : " << nverts << " number of MST edges : " << MST.size() << std::endl;
    //for (int i = 0; i < MST.size(); i++)
    //    std::cout << i << " " << MST[i] << std::endl;
   
    for (int i = 0; i < nEdges(); i++)
    {
        std::cout << "edge : " << i << " " << edgeVertex(i,0) << " " << edgeVertex(i,1) << std::endl;
        //std::cout << "edgeOppositeVertex : " << edgeOppositeVertex(i,0) << " " << edgeOppositeVertex(i,1) << std::endl;
        //std::cout << "edge face : " << edgeFace(i,0) << " " << edgeFace(i, 1) << std::endl;
        for (int j = 0; j < edgeParentID[i].size(); j++)
            std::cout << edgeParentID[i][j] << " ";
        std::cout << std::endl;
    }
    
    return true;   
}
*/
