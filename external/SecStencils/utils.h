#ifndef UTILS_H
#define UTILS_H

#include "types.h"

void 
BuildDiagonalMatrix(
    const VectorX &values,
    SparseMatrixX &A);

void 
BuildSelectorMatrix(
    const size_t colCount,
    const std::vector<int> &rowSelection,
    SparseMatrixX &A);

void
ConvertToVector3(
    const MatrixX& X,
    std::vector<Vector3>& P);

void
ConvertToMatrix(
    const std::vector<Vector3>& P,
    MatrixX& X);

Scalar 
ComputeLinf(const SparseMatrixX& A);

Scalar 
ComputeLi(const VectorX& x);

Scalar
ComputeL2(const SparseMatrixX& A, const VectorX& x);

void
PrintSparse(const SparseMatrixX& A);

// Return index of "query" in "values"
// Return -1 if "query" not found
int 
SearchIndex(
    const std::vector<int>& values, 
    int query);

Scalar
Multiply(const VectorX& x, const SparseMatrixX& A, const VectorX& y);

bool
GetSubdEdgeSign(int index);

#endif // UTILS_H
