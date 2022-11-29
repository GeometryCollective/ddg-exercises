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

// Some sanity checks for the Laplace matrix (same as from poisson-problem)
TEST_F(MeanCurvatureFlowTest, laplaceMatrixSymmetric) {
    bool success = true;
    double eps = 1e-5;
    SparseMatrix<double> A = MCF.geometry->laplaceMatrix();
    size_t V = A.rows();
    for (size_t i = 0; i < V; i++) {
        for (SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
            if (abs(it.value() - A.coeffRef(it.col(), it.row())) > eps) {
                success = false;
                break;
            }
        }
    }
    EXPECT_TRUE(success) << "Laplace matrix is not symmetric";
}

TEST_F(MeanCurvatureFlowTest, laplaceMatrixRowSums) {
    bool success = true;
    double eps = 1e-5;
    SparseMatrix<double> A = MCF.geometry->laplaceMatrix(); // A should equal A^T
    size_t V = A.rows();
    for (size_t i = 0; i < V; i++) {
        double rowSum = 0.;
        for (SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
            rowSum += it.value();
        }
        if (abs(rowSum) > eps) success = false;
    }
    EXPECT_TRUE(success) << "Each row of Laplace matrix does not sum to zero";
}
TEST_F(MeanCurvatureFlowTest, laplaceMatrixDiagonalSums) {
    bool success = true;
    double eps = 1e-5;
    SparseMatrix<double> A = MCF.geometry->laplaceMatrix(); // A should equal A^T
    size_t V = A.rows();
    for (size_t i = 0; i < V; i++) {
        double diag = 0.;
        double sum = 0.;
        for (SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
            if (it.row() != it.col()) {
                sum += it.value();
            } else {
                diag = it.value();
            }
        }
        if (abs(diag + sum) > eps) success = false;
    }
    EXPECT_TRUE(success)
        << "Off-diagonal entries of Laplace matrix do not sum to (negative) the diagonal entry in each row";
}

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