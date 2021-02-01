#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class ScalarPoissonProblem {

  public:
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;
    SparseMatrix<double> A; // Laplace matrix
    SparseMatrix<double> M; // mass matrix
    double totalArea;

    // constructors
    ScalarPoissonProblem() {}
    ScalarPoissonProblem(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);

    Vector<double> solve(const Vector<double>& rho) const;
};