#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <map>

using namespace geometrycentral;
using namespace geometrycentral::surface;

class TreeCotree {

  public:
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;

    std::map<Vertex, Vertex> vertexParent;
    std::map<Face, Face> faceParent;

    std::vector<std::vector<Halfedge>> generators;

    TreeCotree() {}
    TreeCotree(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);

    void buildPrimalSpanningTree();

    bool inPrimalSpanningTree(Halfedge he);

    void buildDualSpanningCoTree();

    bool inDualSpanningCotree(Halfedge he);

    Halfedge sharedHalfedge(Face f, Face g) const;

    void buildGenerators();
};