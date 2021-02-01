#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "solvers.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class HodgeDecomposition {

  public:
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;

    SparseMatrix<double> hodge1;
    SparseMatrix<double> hodge2;
    SparseMatrix<double> d0;
    SparseMatrix<double> d1;
    SparseMatrix<double> hodge1Inv;
    SparseMatrix<double> hodge2Inv;
    SparseMatrix<double> d0T;
    SparseMatrix<double> d1T;
    SparseMatrix<double> A;
    SparseMatrix<double> B;

    HodgeDecomposition() {}
    HodgeDecomposition(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);

    Vector<double> computeExactComponent(const Vector<double>& omega) const;

    Vector<double> computeCoExactComponent(const Vector<double>& omega) const;

    Vector<double> computeHarmonicComponent(const Vector<double>& omega, const Vector<double>& dAlpha,
                                            const Vector<double>& deltaBeta) const;
};