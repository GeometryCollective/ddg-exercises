#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "hodge-decomposition.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class HarmonicBases {

  public:
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;

    HarmonicBases() {}
    HarmonicBases(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);

    Vector<double> buildClosedPrimalOneForm(const std::vector<Halfedge>& generator) const;

    std::vector<Vector<double>> compute(const std::vector<std::vector<Halfedge>>& generators,
                                        const HodgeDecomposition& hodgeDecomposition) const;
};