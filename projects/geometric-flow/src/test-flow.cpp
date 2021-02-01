#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "mean-curvature-flow.h"
#include "modified-mean-curvature-flow.h"
#include "gtest/gtest.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace {

class MeanCurvatureFlowTest : public ::testing::Test {

  protected:
    MeanCurvatureFlow MCF;
    int steps;
    double h;
    std::vector<Vector3> mcf_soln;
    std::vector<Vector3> mmcf_soln;

    void loadSolutionFile() {

        std::vector<Vector3> v;
        std::vector<std::array<int, 3>> f;

        std::string filepath = "../include/test-flow-soln.txt";
        std::ifstream input_file(filepath);
        std::string line;
        std::string X;

        if (input_file.is_open()) {
            while (!input_file.eof()) {
                getline(input_file, line);
                std::istringstream iss(line);
                iss >> X;
                if (X == "v") {
                    double x, y, z;
                    iss >> x >> y >> z;
                    v.push_back({x, y, z});

                } else if (X == "f") {
                    int a, b, c;
                    iss >> a >> b >> c;
                    f.push_back({a, b, c});

                } else if (X == "mcf") {
                    double x, y, z;
                    iss >> x >> y >> z;
                    mcf_soln.push_back({x, y, z});
                } else if (X == "mmcf") {
                    double x, y, z;
                    iss >> x >> y >> z;
                    mmcf_soln.push_back({x, y, z});
                } else if (X == "steps") {
                    int val;
                    iss >> val;
                    steps = val;
                } else if (X == "h") {
                    double val;
                    iss >> val;
                    h = val;
                }
            }
        } else {
            std::cerr << "Oops, could not open input file <" << filepath
                      << "> for unit testing. Make sure filepath is correct." << std::endl;
            std::runtime_error("");
        }

        size_t nVertices = v.size();
        size_t nFaces = f.size();
        Eigen::MatrixXd vMat(nVertices, 3);
        Eigen::MatrixXi fMat(nFaces, 3);
        for (size_t i = 0; i < nVertices; i++) {
            vMat(i, 0) = v[i][0];
            vMat(i, 1) = v[i][1];
            vMat(i, 2) = v[i][2];
        }
        for (size_t i = 0; i < nFaces; i++) {
            fMat(i, 0) = f[i][0] - 1; // Geometry Central takes 0-indexed vertices
            fMat(i, 1) = f[i][1] - 1;
            fMat(i, 2) = f[i][2] - 1;
        }
        std::unique_ptr<ManifoldSurfaceMesh> mesh;
        std::unique_ptr<VertexPositionGeometry> geometry;
        std::tie(mesh, geometry) = makeManifoldSurfaceMeshAndGeometry(vMat, fMat);

        MCF = MeanCurvatureFlow(mesh.release(), geometry.release());
    }

    MeanCurvatureFlowTest() {
        loadSolutionFile();
    }

    virtual ~MeanCurvatureFlowTest() {
        delete MCF.mesh;
        delete MCF.geometry;
    }
};

class ModifiedMeanCurvatureFlowTest : public MeanCurvatureFlowTest {

  protected:
    ModifiedMeanCurvatureFlow mMCF;

    ModifiedMeanCurvatureFlowTest() {
        mMCF = ModifiedMeanCurvatureFlow(MeanCurvatureFlowTest::MCF.mesh, MeanCurvatureFlowTest::MCF.geometry);
    }

    // mesh, geometry already deleted in MCF
};

TEST_F(MeanCurvatureFlowTest, integrate) {

    for (int i = 0; i < steps; i++) {
        MCF.integrate(h);
    }

    MCF.geometry->normalize({0, 0, 0}, false);

    bool success = true;
    for (Vertex v : MCF.mesh->vertices()) {
        if ((MCF.geometry->inputVertexPositions[v] - mcf_soln[v.getIndex()]).norm() > 1e-6) {
            success = false;
            break;
        }
    }
    EXPECT_TRUE(success);
}

TEST_F(ModifiedMeanCurvatureFlowTest, integrate) {

    // don't need to reset mesh positions, since gtest calls constructor with each test
    for (int i = 0; i < steps; i++) {
        mMCF.integrate(h);
    }

    mMCF.geometry->normalize({0, 0, 0}, false);

    bool success = true;
    for (Vertex v : mMCF.mesh->vertices()) {
        if ((mMCF.geometry->inputVertexPositions[v] - mmcf_soln[v.getIndex()]).norm() > 1e-6) {
            success = false;
            break;
        }
    }
    EXPECT_TRUE(success);
}
} // namespace

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}