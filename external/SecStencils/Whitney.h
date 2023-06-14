#ifndef WHITNEY_H
#define WHITNEY_H

#include "Subd.h" 

class Whitney : public Subd
{
public:
    Whitney() { }

    ~Whitney() { }
    
    void BuildS0(SparseMatrixX& A) const;
    void BuildS1(SparseMatrixX& A) const;
    void BuildS2(SparseMatrixX& A) const;

    void GetSubdividedEdges(std::vector< std::vector<int> >& edgeToVert) const;
    void GetSubdividedFaces(std::vector< std::vector<int> >& faceToVert) const;

    bool IsVertRegular(int vert) const;
    bool AreIrregularVertsIsolated() const;

protected:
    int _GetVertVertIndex(int vert) const;
    int _GetEdgeVertIndex(int edge) const;
    int _GetFaceVertIndex(int face) const;

    int _GetEdgeEdgeIndex(int edge, int vertInEdge) const;
    int _GetFaceEdgeIndex(int face, int edgeInFace) const;

    int _GetCentralFaceIndex(int face) const;
    int _GetCornerFaceIndex(int face, int vertInFace) const;
};

#endif // WHITNEY_H
