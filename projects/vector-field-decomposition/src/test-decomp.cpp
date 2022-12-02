#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "harmonic-bases.h"
#include "hodge-decomposition.h"
#include "tree-cotree.h"
#include "gtest/gtest.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace {

class HodgeDecompositionTest : public ::testing::Test {

  protected:
    HodgeDecomposition HD;
    Vector<double> omega;
    Vector<double> dAlpha_soln;
    Vector<double> deltaBeta_soln;
    Vector<double> gamma_soln;

    HodgeDecompositionTest() {

        std::vector<std::array<double, 3>> v;
        std::vector<std::array<int, 3>> f;

        std::vector<double> w;
        std::vector<double> dAlpha;
        std::vector<double> deltaBeta;
        std::vector<double> gamma;

        std::string filepath = "../include/test-decomp-soln.txt";
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

                } else if (X == "omega") {
                    double val;
                    iss >> val;
                    w.push_back(val);

                } else if (X == "dAlpha") {
                    double val;
                    iss >> val;
                    dAlpha.push_back(val);

                } else if (X == "deltaBeta") {
                    double val;
                    iss >> val;
                    deltaBeta.push_back(val);

                } else if (X == "gamma") {
                    double val;
                    iss >> val;
                    gamma.push_back(val);
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
        HD = HodgeDecomposition(mesh.release(), geometry.release());

        size_t nEdges = HD.mesh->nEdges();
        omega = Vector<double>(nEdges);
        dAlpha_soln = Vector<double>(nEdges);
        deltaBeta_soln = Vector<double>(nEdges);
        gamma_soln = Vector<double>(nEdges);
        for (size_t i = 0; i < nEdges; i++) {
            omega[i] = w[i];
            dAlpha_soln[i] = dAlpha[i];
            deltaBeta_soln[i] = deltaBeta[i];
            gamma_soln[i] = gamma[i];
        }
    }

    virtual ~HodgeDecompositionTest() {
        delete HD.mesh;
        delete HD.geometry;
    }
};

class HarmonicBasesTest : public HodgeDecompositionTest {

  protected:
    HarmonicBases HB;

    HarmonicBasesTest() {
        HB = HarmonicBases(HD.mesh, HD.geometry);
    }

    virtual ~HarmonicBasesTest() {
        // mesh, geometry will have already been deleted
    }
};

TEST_F(HodgeDecompositionTest, computeExactComponent) {

    Vector<double> dAlpha = HD.computeExactComponent(this->omega);
    EXPECT_TRUE((dAlpha - dAlpha_soln).norm() < 1e-5);
}

TEST_F(HodgeDecompositionTest, computeCoExactComponent) {

    Vector<double> deltaBeta = HD.computeCoExactComponent(this->omega);
    EXPECT_TRUE((deltaBeta - deltaBeta_soln).norm() < 1e-6);
}

TEST_F(HodgeDecompositionTest, computeHarmonicComponent) {

    Vector<double> gamma = HD.computeHarmonicComponent(this->omega, this->dAlpha_soln, this->deltaBeta_soln);
    EXPECT_TRUE((gamma - gamma_soln).norm() < 1e-6);
}

TEST_F(HarmonicBasesTest, compute) {

    TreeCotree treeCotree = TreeCotree(HB.mesh, HB.geometry);
    treeCotree.buildGenerators();
    std::vector<Vector<double>> bases = HB.compute(treeCotree.generators, HD);
    int N = bases.size();
    Eigen::MatrixXd M(bases[0].size(), N);
    for (int i = 0; i < N; i++) {
        M.col(i) = bases[i];
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(M);
    EXPECT_TRUE(svd.rank() == N) << "Bases are not linearly independent";

    int g = 1 - HB.geometry->eulerCharacteristic() / 2;
    EXPECT_EQ(bases.size(), 2 * g) << "Basis size is not 2g";

    bool closed_success = true;
    bool exact_success = true;
    for (int i = 0; i < N; i++) {
        Vector<double> dGamma = HD.d1 * bases[i];
        Vector<double> deltaGamma = HD.d0T * HD.hodge1 * bases[i];
        if (dGamma.norm() > 1e-4) {
            closed_success = false;
        }
        if (deltaGamma.norm() > 1e-4) {
            exact_success = false;
        }
    }
    EXPECT_TRUE(closed_success) << "Bases are not closed";
    EXPECT_TRUE(exact_success) << "Bases are not exact";
}

} // namespace


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
