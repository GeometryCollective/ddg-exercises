#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "gtest/gtest.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace {

class DiscreteExteriorCalculusTest : public ::testing::Test {

  protected:
    // member variables
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;

    SparseMatrix<double> D0;
    SparseMatrix<double> D1;
    SparseMatrix<double> H0;
    SparseMatrix<double> H1;
    SparseMatrix<double> H2;

    Vector<double> phi;
    Vector<double> omega;
    Vector<double> dPhi_soln;
    Vector<double> dOmega_soln;
    Vector<double> H0_soln;
    Vector<double> H1_soln;
    Vector<double> H2_soln;

    // constructor
    DiscreteExteriorCalculusTest() {

        std::vector<double> d0;
        std::vector<double> d1;
        std::vector<double> h0;
        std::vector<double> h1;
        std::vector<double> h2;
        std::vector<double> p;
        std::vector<double> dPhi;
        std::vector<double> w;
        std::vector<double> dOmega;
        std::vector<std::array<double, 3>> v;
        std::vector<std::array<int, 3>> f;

        std::string filepath = "../include/test-dec-soln.txt";
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

                } else if (X == "hodge0") {
                    double val;
                    iss >> val;
                    h0.push_back(val);

                } else if (X == "phi") {
                    double val;
                    iss >> val;
                    p.push_back(val);

                } else if (X == "hodge1") {
                    double val;
                    iss >> val;
                    h1.push_back(val);

                } else if (X == "dPhi") {
                    double val;
                    iss >> val;
                    dPhi.push_back(val);

                } else if (X == "omega") {
                    double val;
                    iss >> val;
                    w.push_back(val);

                } else if (X == "hodge2") {
                    double val;
                    iss >> val;
                    h2.push_back(val);

                } else if (X == "dOmega") {
                    double val;
                    iss >> val;
                    dOmega.push_back(val);
                }
            }
        } else {
            std::cerr << "Oops, could not open input file <" << filepath
                      << "> for unit testing. Make sure filepath is correct." << std::endl;
            std::runtime_error(""); // errors not enabled
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
        std::tie(this->mesh, this->geometry) = makeManifoldSurfaceMeshAndGeometry(vMat, fMat);

        this->D0 = geometry->buildExteriorDerivative0Form();
        this->D1 = geometry->buildExteriorDerivative1Form();
        this->H0 = geometry->buildHodgeStar0Form();
        this->H1 = geometry->buildHodgeStar1Form();
        this->H2 = geometry->buildHodgeStar2Form();

        this->phi = Vector<double>::Zero(mesh->nVertices());
        this->dPhi_soln = Vector<double>::Zero(mesh->nEdges());
        this->omega = Vector<double>::Zero(mesh->nEdges());
        this->dOmega_soln = Vector<double>::Zero(mesh->nFaces());
        this->H0_soln = Vector<double>::Zero(mesh->nVertices());
        this->H1_soln = Vector<double>::Zero(mesh->nEdges());
        this->H2_soln = Vector<double>::Zero(mesh->nFaces());

        for (size_t i = 0; i < nVertices; i++) {
            phi[i] = p[i];
            H0_soln[i] = h0[i];
        }
        for (size_t i = 0; i < mesh->nEdges(); i++) {
            dPhi_soln[i] = dPhi[i];
            omega[i] = w[i];
            H1_soln[i] = h1[i];
        }
        for (size_t i = 0; i < mesh->nFaces(); i++) {
            dOmega_soln[i] = dOmega[i];
            H2_soln[i] = h2[i];
        }
    }

    virtual ~DiscreteExteriorCalculusTest() {}
};

TEST_F(DiscreteExteriorCalculusTest, buildHodgeStar0Form) {

    std::cerr << "Testing buildHodgeStar0Form()..." << std::endl;
    Vector<double> result = H0 * Vector<double>::Ones(mesh->nVertices());
    EXPECT_TRUE((result - H0_soln).norm() < 1e-6);
}

TEST_F(DiscreteExteriorCalculusTest, buildHodgeStar1Form) {

    std::cerr << "Testing buildHodgeStar1Form()..." << std::endl;
    Vector<double> result = H1 * Vector<double>::Ones(mesh->nEdges());
    EXPECT_TRUE((result - H1_soln).norm() < 1e-4);
}

TEST_F(DiscreteExteriorCalculusTest, buildHodgeStar2Form) {

    std::cerr << "Testing buildHodgeStar2Form()..." << std::endl;
    Vector<double> result = H2 * Vector<double>::Ones(mesh->nFaces());
    EXPECT_TRUE((result - H2_soln).norm() / H2_soln.norm() < 1e-6);
}

TEST_F(DiscreteExteriorCalculusTest, buildExteriorDerivative0Form) {

    std::cerr << "Testing buildExteriorDerivative0Form()..." << std::endl;
    Vector<double> result = D0 * phi;
    // Allow for either choice of edge orientation.
    EXPECT_TRUE((result - dPhi_soln).norm() < 1e-6 || (result + dPhi_soln).norm() < 1e-6);
}

TEST_F(DiscreteExteriorCalculusTest, buildExteriorDerivative1Form) {

    std::cerr << "Testing buildExteriorDerivative1Form()..." << std::endl;
    Vector<double> result = D1 * omega;
    // Allow for either choice of edge orientation.
    EXPECT_TRUE((result - dOmega_soln).norm() < 1e-6 || (result + dOmega_soln).norm() < 1e-6);
}

TEST_F(DiscreteExteriorCalculusTest, d2) {

    std::cerr << "Testing if the discrete exterior derivative applied twice equals 0..." << std::endl;
    SparseMatrix<double> d0d1 = D1 * D0;
    EXPECT_TRUE(d0d1.norm() < 1e-6) << "\t d^2 != 0";
}
} // namespace


int main(int argc, char** argv) {

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}