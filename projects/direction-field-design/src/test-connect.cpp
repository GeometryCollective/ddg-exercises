#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "harmonic-bases.h"
#include "hodge-decomposition.h"
#include "solvers.h"
#include "tree-cotree.h"
#include "trivial-connections.h"
#include "gtest/gtest.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace {

class TrivialConnectionsTest : public ::testing::Test {

  protected:
    TrivialConnections TC;
    Vector<double> singularity;
    Vector<double> deltaBeta_sol;
    Vector<double> gamma_sol;
    Vector<double> phi_sol;
    Vector<double> u;

    TrivialConnectionsTest() {
        intializeSolution("../include/test-connect-soln.txt");
    }

    void intializeSolution(std::string filepath) {

        std::vector<std::array<double, 3>> v;
        std::vector<std::array<int, 3>> f;
        std::vector<double> S;
        std::vector<double> B;
        std::vector<double> G;
        std::vector<double> P;
        std::vector<double> U;

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
                } else if (X == "singularity") {
                    double val;
                    iss >> val;
                    S.push_back(val);
                } else if (X == "deltaBeta") {
                    double val;
                    iss >> val;
                    B.push_back(val);
                } else if (X == "gamma") {
                    double val;
                    iss >> val;
                    G.push_back(val);
                } else if (X == "phi") {
                    double val;
                    iss >> val;
                    P.push_back(val);
                } else if (X == "u") {
                    double val;
                    iss >> val;
                    U.push_back(val);
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
        TC = TrivialConnections(mesh.release(), geometry.release());

        singularity = Vector<double>::Zero(TC.mesh->nVertices());
        for (size_t i = 0; i < S.size(); i++) {
            singularity[i] = S[i];
        }
        if (B.size() > 0) {
            deltaBeta_sol = Vector<double>::Zero(B.size());
            for (size_t i = 0; i < B.size(); i++) {
                deltaBeta_sol[i] = B[i];
            }
        }
        if (G.size() > 0) {
            gamma_sol = Vector<double>::Zero(G.size());
            for (size_t i = 0; i < G.size(); i++) {
                gamma_sol[i] = G[i];
            }
        }
        if (P.size() > 0) {
            phi_sol = Vector<double>::Zero(P.size());
            for (size_t i = 0; i < P.size(); i++) {
                phi_sol[i] = P[i];
            }
        }
        if (U.size() > 0) {
            u = Vector<double>::Zero(U.size());
            for (size_t i = 0; i < U.size(); i++) {
                u[i] = U[i];
            }
        }
    }

    virtual ~TrivialConnectionsTest() {
        delete TC.mesh;
        delete TC.geometry;
    }
};

class TrivialConnectionsSphereTest : public TrivialConnectionsTest {

  protected:
    TrivialConnectionsSphereTest() {
        intializeSolution("../include/test-connect-sphere.txt");
    }

    ~TrivialConnectionsSphereTest() {
        // mesh, geometry will already have been deleted
    }
};

TEST_F(TrivialConnectionsSphereTest, computeConnections) {

    SparseMatrix<double> hodge1 = TC.geometry->buildHodgeStar1Form();
    SparseMatrix<double> hodge1Inv = sparseInverseDiagonal(hodge1);
    SparseMatrix<double> d1 = TC.geometry->buildExteriorDerivative1Form();
    Vector<double> phi = TC.computeConnections(singularity);
    Vector<double> codifferential = d1 * hodge1Inv * phi;
    EXPECT_TRUE(codifferential.norm() < 1e-4) << "Connection form is not co-closed";

    SparseMatrix<double> d0T = TC.geometry->buildExteriorDerivative0Form().transpose();
    Vector<double> differential = d0T * phi;
    EXPECT_TRUE((differential - u).norm() < 1e-4) << "dδβ != u";
}
} // namespace


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}