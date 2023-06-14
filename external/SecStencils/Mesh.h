#ifndef MESH_H
#define MESH_H

#include "types.h"
#include "utils.h"

class Mesh
{
protected:
    bool _isTriangulated;
    bool _isQuadrangulated;
    std::vector<Vector3> _points;
    std::vector< std::vector<int> > _faceToVert; // sorted counter-clockwise
    std::vector< std::vector<int> > _faceToEdge; // ith edge is outgoing from ith vert
    std::vector< std::vector<int> > _edgeToVert; // sorted by increasing index 
    std::vector< std::vector<int> > _edgeToFace; // not sorted
    std::vector< std::vector<int> > _vertToFace; // sorted counter-clockwise
    std::vector< std::vector<int> > _vertToEdge; // ith edge is outgoing from vert in ith face

public:
    Mesh() 
    : _isTriangulated(false)
    , _isQuadrangulated(false)
    { }

    Mesh(const Mesh& mesh) 
    { 
        Copy(mesh); 
    }

    Mesh& operator = (const Mesh& mesh)
    {
        Copy(mesh);
        return *this;
    }

    ~Mesh() { }
    
    inline bool IsTriangulated()   const { return _isTriangulated;   }
    inline bool IsQuadrangulated() const { return _isQuadrangulated; }

    inline size_t GetFaceCount() const { return _faceToVert.size(); }
    inline size_t GetEdgeCount() const { return _edgeToFace.size(); }
    inline size_t GetVertCount() const { return _vertToEdge.size(); }

    inline const std::vector<int>& GetFaceVerts(int face) const { return _faceToVert[face]; }
    inline const std::vector<int>& GetFaceEdges(int face) const { return _faceToEdge[face]; }

    inline const std::vector<int>& GetEdgeVerts(int edge) const { return _edgeToVert[edge]; }
    inline const std::vector<int>& GetEdgeFaces(int edge) const { return _edgeToFace[edge]; }

    inline const std::vector<int>& GetVertEdges(int vert) const { return _vertToEdge[vert]; }
    inline const std::vector<int>& GetVertFaces(int vert) const { return _vertToFace[vert]; }

    inline int GetVertIndexInFace(int face, int vert) const { return SearchIndex(_faceToVert[face], vert); }
    inline int GetEdgeIndexInFace(int face, int edge) const { return SearchIndex(_faceToEdge[face], edge); }

    inline int GetVertIndexInEdge(int edge, int vert) const { return SearchIndex(_edgeToVert[edge], vert); }
    inline int GetFaceIndexInEdge(int edge, int face) const { return SearchIndex(_edgeToFace[edge], face); }

    inline int GetEdgeIndexInVert(int vert, int edge) const { return SearchIndex(_vertToEdge[vert], edge); }
    inline int GetFaceIndexInVert(int vert, int face) const { return SearchIndex(_vertToFace[vert], face); }

    inline bool IsEdgeBoundary(int edge) const { return _edgeToFace[edge].size() == 1; }
    inline bool IsVertBoundary(int vert) const { return _vertToEdge[vert].size() != _vertToFace[vert].size(); }

    inline const Vector3& GetVertPos(int vert) const     { return _points[vert]; }
    inline void SetVertPos(int vert, const Vector3 &pos) { _points[vert] = pos;  }

    inline void GetPos(MatrixX& X) const { ConvertToMatrix (_points, X); }
    inline void SetPos(const MatrixX& X) { ConvertToVector3(X, _points); }

    bool HasBoundary() const;
    
    // Return +1 if face has (edgeVert0,edgeVert1)
    // Return -1 if face has (edgeVert1,edgeVert0)
    int GetEdgeSignInFace(int face, int edgeInFace) const;

    // Return +1 if edge ends at vert
    // Return -1 if edge starts at vert
    int GetVertSignInEdge(int edge, int vertInEdge) const;

    void GetBoundaryVerts(std::vector< std::vector<int> > &verts) const;
    void GetBoundaryFromVert(int vert, std::vector<int>& boundary) const;
    
    void GetBoundaryVerts(std::vector<int> &verts) const;
    void GetBoundaryEdges(std::vector<int> &edges) const;

    void GetNonBoundaryVerts(std::vector<int> &verts) const;
    void GetNonBoundaryEdges(std::vector<int> &edges) const;

    void BuildD(int i, SparseMatrixX& D) const;
    void BuildD0(SparseMatrixX &D) const;
    void BuildD1(SparseMatrixX &D) const;

    ///////////////////
    // Populate mesh //
    ///////////////////

    void Clear();

    void Copy(const Mesh& mesh);

    void Populate(
        const std::vector<Vector3> &points,
        const std::vector<int> &faceVertCount, 
        const std::vector<int> &faceToVert);   // sorted CCW

    void Populate(
        const std::vector<Vector3> &points,
        const std::vector< std::vector<int> > &faceToVert); // sorted CCW

    void Populate(
        const std::vector<Vector3> &points,
        const std::vector< std::vector<int> > &faceToVert,  // sorted CCW
        const std::vector< std::vector<int> > &edgeToVert); // sorted by increasing vert index

    void Triangulate();

protected:
    void _CheckValidity() const;

    void _SortOneRings();

    bool _AreFacesTris() const;

    bool _AreFacesQuads() const;
};

#endif // MESH_H
