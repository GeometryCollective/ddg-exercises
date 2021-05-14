#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "heat-method.h"
#include "gtest/gtest.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace {

class HeatMethodTest : public ::testing::Test {

  protected:
    HeatMethod HM;

    Vector<double> delta;
    Vector<double> phi_soln;
    Vector<double> div_soln;
    FaceData<Vector3> X_soln;

    HeatMethodTest() {

        std::vector<double> d;
        std::vector<double> phi;
        std::vector<double> div;
        std::vector<Vector3> X;

        std::vector<Vector3> v;
        std::vector<std::array<int, 3>> f;

        std::string filepath = "../include/test-heat-soln.txt";
        std::ifstream input_file(filepath);
        std::string line;
        std::string h;

        if (input_file.is_open()) {
            while (!input_file.eof()) {
                getline(input_file, line);
                std::istringstream iss(line);
                iss >> h;
                if (h == "v") {
                    double x, y, z;
                    iss >> x >> y >> z;
                    v.push_back({x, y, z});

                } else if (h == "f") {
                    int a, b, c;
                    iss >> a >> b >> c;
                    f.push_back({a, b, c});

                } else if (h == "delta") {
                    double val;
                    iss >> val;
                    d.push_back(val);

                } else if (h == "X") {
                    double x, y, z;
                    iss >> x >> y >> z;
                    X.push_back(Vector3{x, y, z});

                } else if (h == "div") {
                    double val;
                    iss >> val;
                    div.push_back(-val);

                } else if (h == "phi") {
                    double val;
                    iss >> val;
                    phi.push_back(val);
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
        HM = HeatMethod(mesh.release(), geometry.release());

        delta = Vector<double>::Zero(HM.mesh->nVertices());
        div_soln = Vector<double>::Zero(HM.mesh->nVertices());
        phi_soln = Vector<double>::Zero(HM.mesh->nVertices());
        X_soln = FaceData<Vector3>(*HM.mesh);
        for (size_t i = 0; i < HM.mesh->nVertices(); i++) {
            delta[i] = d[i];
            div_soln[i] = div[i];
            phi_soln[i] = phi[i];
        }
        for (size_t i = 0; i < HM.mesh->nFaces(); i++) {
            X_soln[HM.mesh->face(i)] = X[i];
        }
    }

    virtual ~HeatMethodTest() {
        delete HM.mesh;
        delete HM.geometry;
    }
};

TEST_F(HeatMethodTest, computeVectorField) {

    Eigen::SimplicialLLT<SparseMatrix<double>> llt(HM.F);
    Vector<double> u = llt.solve(delta);
    FaceData<Vector3> X = HM.computeVectorField(u);
    bool success = true;
    for (Face f : HM.mesh->faces()) {
        if ((X[f] - X_soln[f]).norm() > 1e-5) {
            success = false;
            break;
        }
    }
    EXPECT_TRUE(success);
}

TEST_F(HeatMethodTest, computeDivergence) {

    Vector<double> div = HM.computeDivergence(X_soln);
    EXPECT_TRUE((div - div_soln).norm() < 1e-5);
}

TEST_F(HeatMethodTest, compute) {

    Vector<double> phi = HM.compute(delta);
    EXPECT_TRUE((phi - phi_soln).norm() < 1e-5);
}

} // namespace

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
