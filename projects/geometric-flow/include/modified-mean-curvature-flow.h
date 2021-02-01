#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "mean-curvature-flow.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class ModifiedMeanCurvatureFlow : public MeanCurvatureFlow {

  public:
    SparseMatrix<double> A;

    ModifiedMeanCurvatureFlow() {}
    ModifiedMeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);

    SparseMatrix<double> buildFlowOperator(const SparseMatrix<double>& M, double h) const override;
};