#ifndef CATMARK_H
#define CATMARK_H

#include "Subd.h" 

class Catmark : public Subd
{
public:
    Catmark() { }

    void BuildS0(SparseMatrixX& A) const;
    void BuildS1(SparseMatrixX& A) const;
    void BuildS2(SparseMatrixX& A) const;

    void GetSubdividedEdges(std::vector< std::vector<int> >& edgeToVert) const;
    void GetSubdividedFaces(std::vector< std::vector<int> >& faceToVert) const;

    bool IsVertRegular(int vert) const;
    bool AreIrregularVertsIsolated() const;

protected:
    Scalar _GetBeta(int vert) const;
    Scalar _GetGamma(int vert) const;

    int _GetVertVertIndex(int vert) const;
    int _GetEdgeVertIndex(int edge) const;
    int _GetFaceVertIndex(int face) const;

    int _GetEdgeEdgeIndex(int edge, int vertInEdge) const;
    int _GetFaceEdgeIndex(int face, int edgeInFace) const;

    int _GetCentralFaceIndex(int face) const;
    int _GetCornerFaceIndex(int face, int vertInFace) const;

    void _AssembleVertFromFace(int face, TripletInserter out) const;
    void _AssembleVertFromBoundaryVert(int vi, TripletInserter out) const;
    void _AssembleVertFromInteriorVert(int vi, TripletInserter out) const;
    void _AssembleVertFromBoundaryEdge(int edge, TripletInserter out) const;
    void _AssembleVertFromInteriorEdge(int edge, TripletInserter out) const;

    void _AssembleEdgeEvenBoundary(int edge, int vertInEdge, TripletInserter out) const;
    void _AssembleEdgeEvenPartialBoundary(int edge, int vertInEdge, TripletInserter out) const;
    void _AssembleEdgeEvenInterior(int edge, int vertInEdge, TripletInserter out) const;
    void _AssembleEdgeOdd(int face, int edgeInFace, TripletInserter out) const;

    void _AssembleFaceCorner(int face, int vertInFace, TripletInserter out) const;

    void _InsertEdgeVertValue(int row, int col, int vert, int rSign, Scalar val, TripletInserter out) const;
    void _InsertFaceVertValue(int row, int col, int face, int rSign, Scalar val, TripletInserter out) const;

    void _GetFaceWeights(int face, int vertInFace, Scalar& f0, Scalar& f1, Scalar& f2) const;
    void _GetVertWeights(int vert, int edgeInVert, Scalar& f0, Scalar& f1, Scalar& f2) const;
    void _GetBdryWeights(int vert, int edgeInVert, VectorX& sigma, VectorX& xi, VectorX& eta) const;
};

#endif // CATMARK_H
