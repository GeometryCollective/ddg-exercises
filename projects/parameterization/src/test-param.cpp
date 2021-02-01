#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "spectral-conformal-parameterization.h"
#include "gtest/gtest.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace {

class SpectralConformalParameterizationTest : public ::testing::Test {

  protected:
    SpectralConformalParameterization SCP;
    SparseMatrix<std::complex<double>> EC_soln;
    std::vector<Vector2> flattening_soln;

    SpectralConformalParameterizationTest() {

        std::vector<Vector3> v;
        std::vector<std::array<int, 3>> f;
        typedef Eigen::Triplet<std::complex<double>> T;
        std::vector<T> tripletList;

        std::string filepath = "../include/test-param-soln.txt";
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

                } else if (X == "T") {
                    double real, imag, row, col;
                    iss >> real >> imag >> row >> col;
                    tripletList.push_back(T(row, col, std::complex<double>(real, imag)));

                } else if (X == "uv") {
                    double u, v;
                    iss >> u >> v;
                    flattening_soln.push_back({u, v});
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
        size_t N = mesh->nVertices();

        SCP = SpectralConformalParameterization(mesh.release(), geometry.release());

        EC_soln = SparseMatrix<std::complex<double>>(N, N);
        EC_soln.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    virtual ~SpectralConformalParameterizationTest() {
        delete SCP.mesh;
        delete SCP.geometry;
    }
};


TEST_F(SpectralConformalParameterizationTest, buildConformalEnergy) {

    SparseMatrix<std::complex<double>> EC = SCP.buildConformalEnergy();
    EXPECT_TRUE((EC - EC_soln).norm() < 1e-6) << "Conformal energy matrix is incorrect";
}

TEST_F(SpectralConformalParameterizationTest, flatten) {

    VertexData<Vector2> flattening = SCP.flatten();
    // normalize
    Vector2 center = {0.0, 0.0};
    for (Vertex v : SCP.mesh->vertices()) {
        center += flattening[v];
    }
    center /= SCP.mesh->nVertices();
    double radius = 0;
    for (Vertex v : SCP.mesh->vertices()) {
        flattening[v] -= center;
        radius = std::max(radius, flattening[v].norm());
    }
    for (Vertex v : SCP.mesh->vertices()) {
        flattening[v] /= radius;
    }

    // See if there's a rotation between the parameterizations
    Vertex v0 = SCP.mesh->vertex(0);
    std::complex<double> p0 = std::complex<double>(flattening[v0][0], flattening[v0][1]);
    std::complex<double> p1 = std::complex<double>(flattening_soln[0][0], flattening_soln[0][1]);
    std::complex<double> rot = p1 / p0;

    bool success = true;
    for (Vertex v : SCP.mesh->vertices()) {
        size_t i = v.getIndex();
        std::complex<double> p = std::complex<double>(flattening[v][0], flattening[v][1]);
        std::complex<double> s = std::complex<double>(flattening_soln[i][0], flattening_soln[i][1]);
        if (std::norm(rot * p - s) > 1e-6 || std::isnan(std::norm(rot))) {
            success = false;
            break;
        }
    }
    EXPECT_TRUE(success) << "Spectral conformal parameterization is incorrect";
}

} // namespace

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}