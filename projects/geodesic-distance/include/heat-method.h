#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class HeatMethod {

  public:
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;
    SparseMatrix<double> A; // Laplace matrix
    SparseMatrix<double> F; // flow matrix

    HeatMethod() {}
    HeatMethod(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo);
    FaceData<Vector3> computeVectorField(const Vector<double>& u) const;
    Vector<double> computeDivergence(const FaceData<Vector3>& X) const;
    Vector<double> compute(const Vector<double>& delta) const;

  private:
    /*
     * Shift φ such that its minimum value is zero.
     */
    void subtractMinimumDistance(Vector<double>& phi) const {

        // Get min value of the vector <phi>
        int nrows = phi.rows();
        double min = std::numeric_limits<double>::infinity();
        for (int i = 0; i < nrows; i++) {
            min = std::min(min, phi[i]);
        }

        // Re-center
        for (int i = 0; i < nrows; i++) {
            phi[i] -= min;
        }
    }
};