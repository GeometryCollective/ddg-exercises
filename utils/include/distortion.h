#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "colormap.h"
#include <algorithm>
#include <map>

using namespace geometrycentral;
using namespace geometrycentral::surface;

class Distortion {

  public:
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;

    Distortion() {}
    Distortion(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo);

    double computeQuasiConformalError(std::vector<std::array<double, 3>>& colors, VertexData<Vector2>& flattening);
    double computeAreaScaling(std::vector<std::array<double, 3>>& colors, VertexData<Vector2>& flattening);

  private:
    double computeQuasiConformalErrorPerFace(Face f, VertexData<Vector2>& flattening);
    double computeAreaScalingPerFace(Face f, VertexData<Vector2>& flattening);
};