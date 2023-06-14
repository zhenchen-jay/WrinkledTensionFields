#include "Catmark.h"

#include <cassert>
#include <iostream>

bool
Catmark::IsVertRegular(int vert) const
{
    assert(_meshPtr);
    if (_meshPtr->IsVertBoundary(vert)) return true;
    return (_meshPtr->GetVertEdges(vert).size() == 4);
}

bool 
Catmark::AreIrregularVertsIsolated() const
{
    assert(_meshPtr);
    for (int edge = 0; edge < _meshPtr->GetEdgeCount(); ++edge)
    {
        const std::vector<int>& eVerts = _meshPtr->GetEdgeVerts(edge);
        if (IsVertRegular(eVerts[0])) continue;
        if (IsVertRegular(eVerts[1])) continue;
        return false;
    }
    return true;
}

int 
Catmark::_GetVertVertIndex(int vert) const
{
    return vert;
}

int 
Catmark::_GetEdgeVertIndex(int edge) const
{
    assert(_meshPtr);
    return _meshPtr->GetVertCount() + edge;    
}

int 
Catmark::_GetFaceVertIndex(int face) const
{
    assert(_meshPtr);
    return _meshPtr->GetVertCount() + _meshPtr->GetEdgeCount() + face;    
}

int 
Catmark::_GetEdgeEdgeIndex(int edge, int vertInEdge) const
{
    return 2*edge + vertInEdge;
}

int 
Catmark::_GetFaceEdgeIndex(int face, int edgeInFace) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsQuadrangulated());
    return 2*_meshPtr->GetEdgeCount() + 4*face + edgeInFace;
}

int 
Catmark::_GetCentralFaceIndex(int face) const
{
    assert(false);
    return -1;
}

int 
Catmark::_GetCornerFaceIndex(int face, int vertInFace) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsQuadrangulated());
    return 4*face + vertInFace;    
}

void 
Catmark::GetSubdividedEdges(std::vector< std::vector<int> >& edgeToVert) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsQuadrangulated());

    int E = _meshPtr->GetEdgeCount();
    int F = _meshPtr->GetFaceCount();
    edgeToVert.resize(2*E + 4*F);

    for (int edge = 0; edge < E; ++edge)
    {
        const std::vector<int>& eVerts = _meshPtr->GetEdgeVerts(edge);
        for (int i = 0; i < 2; ++i)
        {
            int v0 = _GetVertVertIndex(eVerts[i]);
            int v1 = _GetEdgeVertIndex(edge);
            if (v0 > v1) std::swap(v0, v1);

            int index = _GetEdgeEdgeIndex(edge, i);
            edgeToVert[index].push_back(v0);
            edgeToVert[index].push_back(v1);
        }
    }

    for (int face = 0; face < F; ++face)
    {
        const std::vector<int>& fEdges = _meshPtr->GetFaceEdges(face);
        for (int i = 0; i < 4; ++i)
        {
            int v0 = _GetEdgeVertIndex(fEdges[i]);
            int v1 = _GetFaceVertIndex(face);
            if (v0 > v1) std::swap(v0, v1);

            int index = _GetFaceEdgeIndex(face, i);
            edgeToVert[index].push_back(v0);
            edgeToVert[index].push_back(v1);
        }
    }
}

void 
Catmark::GetSubdividedFaces(std::vector< std::vector<int> >& faceToVert) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsQuadrangulated());

    int V = _meshPtr->GetVertCount();
    int F = _meshPtr->GetFaceCount();
    faceToVert.resize(4*F);

    for (int face = 0; face < _meshPtr->GetFaceCount(); ++face)
    {
        const std::vector<int>& fVerts = _meshPtr->GetFaceVerts(face);
        const std::vector<int>& fEdges = _meshPtr->GetFaceEdges(face);
        for (int j = 0; j < 4; ++j)
        {
            int index = _GetCornerFaceIndex(face, j);
            faceToVert[index].push_back(_GetVertVertIndex(fVerts[j]));
            faceToVert[index].push_back(_GetEdgeVertIndex(fEdges[j]));
            faceToVert[index].push_back(_GetFaceVertIndex(face));
            faceToVert[index].push_back(_GetEdgeVertIndex(fEdges[(j+3)%4]));
        }
    }    
}

Scalar
Catmark::_GetBeta(int vert) const
{
    // Page 130 in [Wang thesis]
    Scalar k = _meshPtr->GetVertFaces(vert).size();
    return 1.5/k;
}

Scalar 
Catmark::_GetGamma(int vert) const
{
    // Page 130 in [Wang thesis]
    Scalar k = _meshPtr->GetVertFaces(vert).size();
    return 0.25/k;
}

void 
Catmark::BuildS0(SparseMatrixX& A) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsQuadrangulated());

    std::vector<TripletX> triplets;
    int V = _meshPtr->GetVertCount();
    int E = _meshPtr->GetEdgeCount();
    int F = _meshPtr->GetFaceCount();

    for (int vi = 0; vi < V; ++vi)
    {
        if (_meshPtr->IsVertBoundary(vi))
            _AssembleVertFromBoundaryVert(vi, std::back_inserter(triplets));
        else
            _AssembleVertFromInteriorVert(vi, std::back_inserter(triplets));
    }

    for (int edge = 0; edge < E; ++edge)
    {
        if (_meshPtr->IsEdgeBoundary(edge))
            _AssembleVertFromBoundaryEdge(edge, std::back_inserter(triplets));
        else 
            _AssembleVertFromInteriorEdge(edge, std::back_inserter(triplets));
    }

    for (int face = 0; face < F; ++face)
    {
        _AssembleVertFromFace(face, std::back_inserter(triplets));
    }

    A.resize(V+E+F, V);
    A.setFromTriplets(triplets.begin(), triplets.end());
}

void
Catmark::_AssembleVertFromFace(int face, TripletInserter out) const
{
    // Fig. C2 (top-left) in [Wang thesis]
    int row = _GetFaceVertIndex(face);
    const std::vector<int>& fVerts = _meshPtr->GetFaceVerts(face);
    for (int i = 0; i < fVerts.size(); ++i)
    {
        *out++ = TripletX(row, fVerts[i], 0.25);
    }
}

void 
Catmark::_AssembleVertFromBoundaryVert(int vi, TripletInserter out) const
{
    // Eq. 9 in [DeRose et al. 1998]
    std::vector<int> boundary(2);
    boundary[0] = _meshPtr->GetVertEdges(vi).front();
    boundary[1] = _meshPtr->GetVertEdges(vi).back();
    
    int row = _GetVertVertIndex(vi);
    for (int j = 0; j < boundary.size(); ++j)
    {
        int edge = boundary[j];
        assert(_meshPtr->IsEdgeBoundary(edge));
        int viInEdge = _meshPtr->GetVertIndexInEdge(edge, vi);
        int vj = _meshPtr->GetEdgeVerts(edge)[(viInEdge+1)%2];
        *out++ = TripletX(row, vj, 0.125);
    }
    *out++ = TripletX(row, vi, 0.75);
}

void 
Catmark::_AssembleVertFromInteriorVert(int vi, TripletInserter out) const
{
    // Fig. C2 (bot-left) in [Wang thesis]
    int row = _GetVertVertIndex(vi);

    const std::vector<int>& vFaces = _meshPtr->GetVertFaces(vi);
    const std::vector<int>& vEdges = _meshPtr->GetVertEdges(vi);
    
    Scalar k = vFaces.size();
    Scalar beta = _GetBeta(vi);
    Scalar gamma = _GetGamma(vi);

    *out++ = TripletX(row, vi, 1.-beta-gamma);

    for (int i = 0; i < vEdges.size(); ++i)
    {
        int edge = vEdges[i];
        int viInEdge = _meshPtr->GetVertIndexInEdge(edge, vi);
        int vj = _meshPtr->GetEdgeVerts(edge)[(viInEdge+1)%2];
        *out++ = TripletX(row, vj, beta/k);
    }

    for (int i = 0; i < vFaces.size(); ++i)
    {
        int face = vFaces[i];
        int viInFace = _meshPtr->GetVertIndexInFace(face, vi);
        int vj = _meshPtr->GetFaceVerts(face)[(viInFace+2)%4];
        *out++ = TripletX(row, vj, gamma/k);
    }
}

void 
Catmark::_AssembleVertFromBoundaryEdge(int edge, TripletInserter out) const
{
    // Eq. 8 in [DeRose et al. 1998]
    int row = _GetEdgeVertIndex(edge);
    const std::vector<int>& eVerts = _meshPtr->GetEdgeVerts(edge);
    for (int i = 0; i < eVerts.size(); ++i)
    {
        *out++ = TripletX(row, eVerts[i], 0.5);
    }
}

void 
Catmark::_AssembleVertFromInteriorEdge(int edge, TripletInserter out) const
{
    // Fig. C2 (top-mid) in [Wang thesis]
    int row = _GetEdgeVertIndex(edge);

    const std::vector<int>& eVerts = _meshPtr->GetEdgeVerts(edge);
    for (int i = 0; i < eVerts.size(); ++i)
    {
        *out++ = TripletX(row, eVerts[i], 0.25);
    }

    const std::vector<int>& eFaces = _meshPtr->GetEdgeFaces(edge);
    for (int i = 0; i < eFaces.size(); ++i)
    {
        const std::vector<int>& fVerts = _meshPtr->GetFaceVerts(eFaces[i]);
        for (int j = 0; j < fVerts.size(); ++j)
        {
            *out++ = TripletX(row, fVerts[j], 0.0625);
        }
    }
}

void 
Catmark::BuildS1(SparseMatrixX& A) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsQuadrangulated());

    std::vector<TripletX> triplets;
    int E = _meshPtr->GetEdgeCount();
    int F = _meshPtr->GetFaceCount();

    for (int edge = 0; edge < E; ++edge)
    {
        if (_meshPtr->IsEdgeBoundary(edge))
        {
            _AssembleEdgeEvenBoundary(edge, 0, std::back_inserter(triplets));
            _AssembleEdgeEvenBoundary(edge, 1, std::back_inserter(triplets));           
        }
        else
        {
            const std::vector<int>& eVerts = _meshPtr->GetEdgeVerts(edge);
            for (int i = 0; i < eVerts.size(); ++i)
            {
                if (_meshPtr->IsVertBoundary(eVerts[i]))
                    _AssembleEdgeEvenPartialBoundary(edge, i, std::back_inserter(triplets));
                else
                    _AssembleEdgeEvenInterior(edge, i, std::back_inserter(triplets));
            }
        }
    }

    for (int face = 0; face < F; ++face)
    {
        _AssembleEdgeOdd(face, 0, std::back_inserter(triplets));
        _AssembleEdgeOdd(face, 1, std::back_inserter(triplets));
        _AssembleEdgeOdd(face, 2, std::back_inserter(triplets));
        _AssembleEdgeOdd(face, 3, std::back_inserter(triplets));
    }

    A.resize(2*E + 4*F, E);
    A.setFromTriplets(triplets.begin(), triplets.end());
}

void
Catmark::_AssembleEdgeEvenPartialBoundary(int edge, int vertInEdge, TripletInserter out) const
{
    //
    // Not covered in [Wang thesis]
    //
    int row = _GetEdgeEdgeIndex(edge, vertInEdge);
    int vert = _meshPtr->GetEdgeVerts(edge)[vertInEdge];
    int rSign = (_GetVertVertIndex(vert) < _GetEdgeVertIndex(edge))? 1 : -1;

    const std::vector<int>& vEdges = _meshPtr->GetVertEdges(vert);
    int edgeCount = vEdges.size();

    int edgeInVert = _meshPtr->GetEdgeIndexInVert(vert, edge);
    assert(edgeInVert != edgeCount-1);
    assert(edgeInVert != 0);

    const std::vector<int>& vFaces = _meshPtr->GetVertFaces(vert);
    int faceCount = vFaces.size();
    assert(faceCount > 1);

    std::vector< std::pair<int, Scalar> > eValues;
    std::vector<TripletX> fValues;

    if (faceCount == 2)
    {
        Scalar f0, f1, f2;
        _GetVertWeights(vert, 1, f0, f1, f2);
        eValues.push_back(std::make_pair(vEdges[1], 0.25 + f0 - f1));
        fValues.push_back(TripletX(vFaces[0], 1,  f1));
        fValues.push_back(TripletX(vFaces[1], 2, -f1));
    }
    else if (edgeInVert == 1)
    {
        Scalar f0, f1, f2;
        _GetVertWeights(vert, 0, f0, f1, f2);

        eValues.push_back(std::make_pair(vEdges[0], 0.1875 - f0));
        eValues.push_back(std::make_pair(vEdges[1], 0.25   + f0 - f1));
        eValues.push_back(std::make_pair(vEdges[2], 0.0625 + f1 - f2));
        eValues.push_back(std::make_pair(vEdges[edgeCount-1], f2 - 0.125));

        fValues.push_back(TripletX(vFaces[0], 1, 0.25   - f0));
        fValues.push_back(TripletX(vFaces[0], 2, 0.1875 - f0));

        fValues.push_back(TripletX(vFaces[1], 1, 0.0625 - f1));
        fValues.push_back(TripletX(vFaces[1], 2,        - f1));
        
        for (int i = 2; i < faceCount; ++i)
        {
            fValues.push_back(TripletX(vFaces[i], 1, -f2));
            fValues.push_back(TripletX(vFaces[i], 2, -f2));
        }
    }
    else if (edgeInVert == edgeCount-2)
    {
        // Symmetric case of edgeInVert==1
        Scalar f0, f1, f2;
        _GetVertWeights(vert, edgeCount-2, f0, f1, f2);

        eValues.push_back(std::make_pair(vEdges[edgeCount-1], 0.1875 - f0));
        eValues.push_back(std::make_pair(vEdges[edgeCount-2], 0.25   + f0 - f1));
        eValues.push_back(std::make_pair(vEdges[edgeCount-3], 0.0625 + f1 - f2));
        eValues.push_back(std::make_pair(vEdges[0], f2 - 0.125));

        fValues.push_back(TripletX(vFaces[faceCount-1], 1, f0 - 0.1875));
        fValues.push_back(TripletX(vFaces[faceCount-1], 2, f0 - 0.25  ));
        
        fValues.push_back(TripletX(vFaces[faceCount-2], 1, f1         ));
        fValues.push_back(TripletX(vFaces[faceCount-2], 2, f1 - 0.0625));
        
        for (int i = 0; i < faceCount-2; ++i)
        {
            fValues.push_back(TripletX(vFaces[i], 1, f2));
            fValues.push_back(TripletX(vFaces[i], 2, f2));
        }
    }
    else
    {
        assert(faceCount > 3);

        VectorX sigma, xi, eta;
        _GetBdryWeights(vert, edgeInVert, sigma, xi, eta);

        eValues.push_back(std::make_pair(vEdges.front(), sigma[0]));
        eValues.push_back(std::make_pair(vEdges.back (), sigma[5]));

        eValues.push_back(std::make_pair(vEdges[edgeInVert-1], sigma[1]));
        eValues.push_back(std::make_pair(vEdges[edgeInVert  ], sigma[2]));
        eValues.push_back(std::make_pair(vEdges[edgeInVert+1], sigma[3]));

        if ((faceCount > 4) && (edgeInVert < edgeCount-3))
        {
            // Additional entry for internal edge of boundary vert
            eValues.push_back(std::make_pair(vEdges[edgeInVert+2], sigma[4]));
        }

        for (int i = 0; i < edgeInVert-1; ++i)
        {
            fValues.push_back(TripletX(vFaces[i], 1, xi [0]));
            fValues.push_back(TripletX(vFaces[i], 2, eta[0]));            
        }

        for (int i = -1; i < 2; ++i)
        {
            fValues.push_back(TripletX(vFaces[edgeInVert+i], 1, xi [2+i]));
            fValues.push_back(TripletX(vFaces[edgeInVert+i], 2, eta[2+i]));
        }

        for (int i = edgeInVert+2; i < faceCount; ++i)
        {
            fValues.push_back(TripletX(vFaces[i], 1, xi [4]));
            fValues.push_back(TripletX(vFaces[i], 2, eta[4]));            
        }
    }

    for (size_t i = 0; i < eValues.size(); ++i)
    {
        _InsertEdgeVertValue(row, eValues[i].first, vert, rSign, eValues[i].second, out);
    }

    for (size_t i = 0; i < fValues.size(); ++i)
    {
        int face = fValues[i].row();
        int vertInFace = _meshPtr->GetVertIndexInFace(face, vert);
        int colInFace = (vertInFace + fValues[i].col()) % 4;
        int col = _meshPtr->GetFaceEdges(face)[colInFace];
        _InsertFaceVertValue(row, col, face, rSign, fValues[i].value(), out);
    }
}

void
Catmark::_AssembleEdgeEvenInterior(int edge, int vertInEdge, TripletInserter out) const
{
    int vert = _meshPtr->GetEdgeVerts(edge)[vertInEdge];
    int edgeInVert = _meshPtr->GetEdgeIndexInVert(vert, edge);

    int row = _GetEdgeEdgeIndex(edge, vertInEdge);
    int rSign = (_GetVertVertIndex(vert) < _GetEdgeVertIndex(edge))? 1 : -1;

    const std::vector<int>& vEdges = _meshPtr->GetVertEdges(vert);
    const std::vector<int>& vFaces = _meshPtr->GetVertFaces(vert);
    
    const int k = vFaces.size();
    assert(k > 2);

    // Page 118 in [Wang thesis]
    const Scalar beta = _GetBeta(vert);
    const Scalar gamma = _GetGamma(vert);

    // Page 118 in [Wang thesis]
    const Scalar alpha0 = 0.5 + gamma;
    const Scalar alpha1 = 0.125 + gamma;
    const Scalar alpha2 = gamma;

    // Page 118 in [Wang thesis]
    VectorX xi = VectorX::Zero(k);
    xi[0] = - 0.0625 - 0.5*gamma/k + 0.125*alpha0;
    xi[1] = xi[0] - 0.0625 + 0.25*alpha1;
    for (int i = 2; i < k-1; ++i) xi[i] = xi[i-1] + 0.25*alpha2;
    xi[k-1] = xi[k-2] + 0.25*alpha1;

    // Page 118 in [Wang thesis]
    VectorX sigma = VectorX::Zero(k);
    sigma[0] = 0.375 - beta/k + 2.*xi[0];
    sigma[1] = 0.0625 - beta/k + xi[1] + xi[k-1];
    for (int i = 2; i < k-1; ++i) sigma[i] = - beta/k + xi[i] + xi[k-i];
    sigma[k-1] = sigma[1];

    for (int i = 0; i < vEdges.size(); ++i)
    {
        _InsertEdgeVertValue(row, vEdges[(edgeInVert+i)%k], vert, rSign, sigma[i], out);
    }

    for (int i = 0; i < vFaces.size(); ++i)
    {
        int face = vFaces[(edgeInVert+i)%k];
        int vertInFace = _meshPtr->GetVertIndexInFace(face, vert);
        const std::vector<int>& fEdges = _meshPtr->GetFaceEdges(face);
        _InsertFaceVertValue(row, fEdges[(vertInFace+1)%4], face, rSign,  xi[    i], out);
        _InsertFaceVertValue(row, fEdges[(vertInFace+2)%4], face, rSign, -xi[k-1-i], out);
    }
}

void
Catmark::_AssembleEdgeEvenBoundary(int edge, int vertInEdge, TripletInserter out) const
{
    //
    // Not covered in [Wang thesis], but same rule as in Loop
    //
    int row = _GetEdgeEdgeIndex(edge, vertInEdge);
    int vert = _meshPtr->GetEdgeVerts(edge)[vertInEdge];
    int rSign = (_GetVertVertIndex(vert) < _GetEdgeVertIndex(edge))? 1 : -1;

    int nEdge = _meshPtr->GetVertEdges(vert).front();
    int vertInNedge = _meshPtr->GetVertIndexInEdge(nEdge, vert);
    int nSign = _meshPtr->GetVertSignInEdge(nEdge, vertInNedge);

    int pEdge = _meshPtr->GetVertEdges(vert).back();
    int vertInPedge = _meshPtr->GetVertIndexInEdge(pEdge, vert);
    int pSign = _meshPtr->GetVertSignInEdge(pEdge, vertInPedge);

    assert(edge == nEdge || edge == pEdge);

    if (edge == nEdge)
    {
        *out++ = TripletX(row, nEdge, (nSign == rSign)? -0.375 :  0.375);
        *out++ = TripletX(row, pEdge, (pSign == rSign)?  0.125 : -0.125);
    }
    else
    {
        *out++ = TripletX(row, pEdge, (pSign == rSign)? -0.375 :  0.375);
        *out++ = TripletX(row, nEdge, (nSign == rSign)?  0.125 : -0.125);
    }
}

void
Catmark::_AssembleEdgeOdd(int face, int edgeInFace, TripletInserter out) const
{
    int row = _GetFaceEdgeIndex(face, edgeInFace);
    int edge = _meshPtr->GetFaceEdges(face)[edgeInFace];
    int rSign = (_GetEdgeVertIndex(edge) < _GetFaceVertIndex(face))? 1 : -1;

    // next edge in face
    int nEdge = _meshPtr->GetFaceEdges(face)[(edgeInFace+1)%4];
    int nSign = _meshPtr->GetEdgeSignInFace(face, (edgeInFace+1)%4);
    
    // prev edge in face
    int pEdge = _meshPtr->GetFaceEdges(face)[(edgeInFace+3)%4];
    int pSign = _meshPtr->GetEdgeSignInFace(face, (edgeInFace+3)%4);

    if (_meshPtr->IsEdgeBoundary(edge))
    {
        // Boundary case not covered in [Wang thesis]
        *out++ = TripletX(row, nEdge, (nSign == rSign)?  0.25 : -0.25);
        *out++ = TripletX(row, pEdge, (pSign == rSign)? -0.25 :  0.25);
    }
    else
    {
        //
        // Fig. C2 (top-right) in [Wang thesis]
        //
        int faceInEdge = _meshPtr->GetFaceIndexInEdge(edge, face);
        int oFace = _meshPtr->GetEdgeFaces(edge)[(faceInEdge+1)%2];
        int edgeInOface = _meshPtr->GetEdgeIndexInFace(oFace, edge);

        int nOppEdge = _meshPtr->GetFaceEdges(oFace)[(edgeInOface+1)%4];
        int nOppSign = _meshPtr->GetEdgeSignInFace(oFace, (edgeInOface+1)%4);

        int pOppEdge = _meshPtr->GetFaceEdges(oFace)[(edgeInOface+3)%4];
        int pOppSign = _meshPtr->GetEdgeSignInFace(oFace, (edgeInOface+3)%4);

        *out++ = TripletX(row, nEdge, (nSign == rSign)?  0.1875 : -0.1875);
        *out++ = TripletX(row, pEdge, (pSign == rSign)? -0.1875 :  0.1875);
    
        *out++ = TripletX(row, nOppEdge, (nOppSign == rSign)? -0.0625 :  0.0625);
        *out++ = TripletX(row, pOppEdge, (pOppSign == rSign)?  0.0625 : -0.0625);
    }
}

void 
Catmark::BuildS2(SparseMatrixX& A) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsQuadrangulated());

    std::vector<TripletX> triplets;
    int F = _meshPtr->GetFaceCount();

    for (int face = 0; face < F; ++face)
    {
        _AssembleFaceCorner(face, 0, std::back_inserter(triplets));
        _AssembleFaceCorner(face, 1, std::back_inserter(triplets));
        _AssembleFaceCorner(face, 2, std::back_inserter(triplets));
        _AssembleFaceCorner(face, 3, std::back_inserter(triplets));
    }

    A.resize(4*F, F);
    A.setFromTriplets(triplets.begin(), triplets.end());
}

void
Catmark::_AssembleFaceCorner(int face, int vertInFace, TripletInserter out) const
{
    //
    // Not covered in [Wang thesis].
    //
    int row = _GetCornerFaceIndex(face, vertInFace);
    int vert = _meshPtr->GetFaceVerts(face)[vertInFace];
    int faceInVert = _meshPtr->GetFaceIndexInVert(vert, face);

    const std::vector<int>& vFaces = _meshPtr->GetVertFaces(vert);
    int k = vFaces.size();

    Scalar f0, f1, f2;
    _GetFaceWeights(face, vertInFace, f0, f1, f2);

    if (_meshPtr->IsVertBoundary(vert))
    {
        if (faceInVert < k-1)
        {
            *out++ = TripletX(row, vFaces[faceInVert+1], f1);
        }
        for (int i = faceInVert+2; i < k; ++i) 
        {
            *out++ = TripletX(row, vFaces[i], f2);
        }
        for (int i = faceInVert-2; i >= 0; --i) 
        {
            *out++ = TripletX(row, vFaces[i], f2);
        }
        if (faceInVert > 0)
        {
            *out++ = TripletX(row, vFaces[faceInVert-1], f1);
        }
        *out++ = TripletX(row, face, f0);
    }
    else
    {
        assert(k > 2);
        *out++ = TripletX(row, vFaces[(faceInVert+1)%k], f1);
        for (int i = 2; i < k-1; ++i) 
        {
            *out++ = TripletX(row, vFaces[(faceInVert+i)%k], f2);
        }
        *out++ = TripletX(row, vFaces[(faceInVert+k-1)%k], f1);
        *out++ = TripletX(row, face, f0);
    }
}

void
Catmark::_GetVertWeights(int vert, int edgeInVert, Scalar& f0, Scalar& f1, Scalar& f2) const
{
    //
    // Subdivision weights for corner face to the left of (vert, edgeInVert).
    //
    int faceInVert = edgeInVert;
    int face = _meshPtr->GetVertFaces(vert)[faceInVert];
    int vertInFace = _meshPtr->GetVertIndexInFace(face, vert);
    _GetFaceWeights(face, vertInFace, f0, f1, f2);
}

void
Catmark::_GetFaceWeights(int face, int vertInFace, Scalar& f0, Scalar& f1, Scalar& f2) const
{
    //
    // Subdivision weights for corner face at (face, vertInFace).
    //
    int vert = _meshPtr->GetFaceVerts(face)[vertInFace];
    const std::vector<int>& vFaces = _meshPtr->GetVertFaces(vert);

    const Scalar gamma = _GetGamma(vert);
    const Scalar alpha1 = 0.125 + gamma;
    const Scalar alpha2 = gamma;

    Scalar sum = 0.;
    int k = vFaces.size();
    if (_meshPtr->IsVertBoundary(vert))
    {
        int faceInVert = _meshPtr->GetFaceIndexInVert(vert, face);
        if (faceInVert > 0)   sum += alpha1;
        if (faceInVert > 1)   sum += (faceInVert - 1)*alpha2;
        if (faceInVert < k-2) sum += (k - 2 - faceInVert)*alpha2;
        if (faceInVert < k-1) sum += alpha1;
    }
    else
    {
        sum = 2*alpha1 + (k-3)*alpha2;
    }

    // Must sum to 0.25.
    f0 = 0.25*(1. - sum);
    f1 = 0.25*alpha1;
    f2 = 0.25*alpha2;
}

void
Catmark::_InsertEdgeVertValue(int row, int col, int vert, int rSign, Scalar val, TripletInserter out) const
{
    // Handy function that sets the sign of val for an edge col incident to vert.
    int vertInCol = _meshPtr->GetVertIndexInEdge(col, vert);
    int sign = _meshPtr->GetVertSignInEdge(col, vertInCol);
    *out++ = TripletX(row, col, (sign == rSign)? -val : val);
}

void
Catmark::_InsertFaceVertValue(int row, int col, int face, int rSign, Scalar val, TripletInserter out) const
{
    // Handy function that sets the sign of val for an edge col incident to face.
    int colInFace = _meshPtr->GetEdgeIndexInFace(face, col);
    int sign = _meshPtr->GetEdgeSignInFace(face, colInFace);
    *out++ = TripletX(row, col, (sign == rSign)? val : -val);
}

void
Catmark::_GetBdryWeights(int vert, int edgeInVert, VectorX& sigma, VectorX& xi, VectorX& eta) const
{
    sigma = VectorX::Zero(6); // [           0, edgeInVert-1, edgeInVert, edgeInVert+1, edgeInVert+2, edgeCount-1]
    xi    = VectorX::Zero(5); // [edgeInVert-2, edgeInVert-1, edgeInVert, edgeInVert+1, edgeInVert+2]
    eta   = VectorX::Zero(5); // [edgeInVert-2, edgeInVert-1, edgeInVert, edgeInVert+1, edgeInVert+2]

    int edgeCount = _meshPtr->GetVertEdges(vert).size();
    int faceCount = _meshPtr->GetVertFaces(vert).size();

    assert(faceCount+1 == edgeCount);
    assert(edgeInVert < edgeCount-2);
    assert(edgeInVert > 1);

    if (edgeInVert == edgeCount-3) // base case
    {
        assert(faceCount > 3);

        Scalar g0, g1, g2;
        _GetVertWeights(vert, edgeCount-2, g0, g1, g2);

        Scalar f0, f1, f2;
        _GetVertWeights(vert, edgeInVert, f0, f1, f2);

        sigma[0] = -0.125  + f2 + g2;
        sigma[1] =  0.0625 + f1 - f2;
        sigma[2] =  0.25   + f0 - f1 + g1 - g2;
        sigma[3] =  0.0625 + f1 - f0 + g0 - g1;
        sigma[4] =  0.;
        sigma[5] =  0.125  - f1 - g0;

        xi[0] = f2 + g2;
        xi[1] = f1 + g2;
        xi[2] = f0 + g1 - 0.1875;
        xi[3] = f1 + g0 - 0.25;
        xi[4] = 0.;
    }
    else // recursive case
    {
        assert(faceCount > 4);

        Scalar f0, f1, f2;
        _GetVertWeights(vert, edgeInVert, f0, f1, f2);

        VectorX pSigma, pXi, pEta;
        _GetBdryWeights(vert, edgeInVert+1, pSigma, pXi, pEta);

        sigma[0] =  f2 + pSigma[0];
        sigma[1] =  f1 - f2             + 0.0625;
        sigma[2] =  f0 - f1 + pSigma[1] + 0.1875;
        sigma[3] =  f1 - f0 + pSigma[2] - 0.1875;
        sigma[4] =  f2 - f1 + pSigma[3] - 0.0625;
        sigma[5] = -f2 + pSigma[5];

        xi[0] = f2 + pXi[0];
        xi[1] = f1 + pXi[0];
        xi[2] = f0 + pXi[1] - 0.1875;
        xi[3] = f1 + pXi[2] - 0.0625;
        xi[4] = f2 + pXi[3];
    }

    eta[0] = xi[0];
    eta[1] = xi[1] - 0.0625;
    eta[2] = xi[2] - 0.0625;
    eta[3] = xi[3];
    eta[4] = xi[4];
}
