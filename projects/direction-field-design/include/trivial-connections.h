#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "harmonic-bases.h"
#include "hodge-decomposition.h"
#include "tree-cotree.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class TrivialConnections {

  public:
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;
    std::vector<Vector<double>> bases;
    SparseMatrix<double> P;
    SparseMatrix<double> A;
    SparseMatrix<double> hodge1;
    SparseMatrix<double> d0;

    TrivialConnections() {}
    TrivialConnections(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);

    SparseMatrix<double> buildPeriodMatrix() const;
    bool satsifyGaussBonnet(const Vector<double>& singularity) const;
    Vector<double> computeCoExactComponent(const Vector<double>& singularity) const;
    double transportNoRotation(Halfedge he, double alphaI) const;
    Vector<double> computeHarmonicComponent(const Vector<double>& deltaBeta) const;
    Vector<double> computeConnections(const Vector<double>& singularity) const;
};