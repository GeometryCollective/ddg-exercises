#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "solvers.h"

#include <complex>

using namespace geometrycentral;
using namespace geometrycentral::surface;

class SpectralConformalParameterization {

  public:
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;

    SpectralConformalParameterization() {}
    SpectralConformalParameterization(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);

    VertexData<Vector2> flatten() const;

    SparseMatrix<std::complex<double>> buildConformalEnergy() const;
};