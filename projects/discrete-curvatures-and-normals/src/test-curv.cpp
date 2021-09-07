#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "gtest/gtest.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace {

class DiscreteCurvaturesAndNormalsTest : public ::testing::Test {

  protected:
    // Member variables
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;

    std::vector<double> angles;
    std::vector<double> dihedralAngles;
    std::vector<Vector3> angleWeightedNormals;
    std::vector<Vector3> sphereInscribedNormals;
    std::vector<Vector3> areaWeightedNormals;
    std::vector<Vector3> gaussCurvatureNormals;
    std::vector<Vector3> meanCurvatureNormals;
    std::vector<double> angleDefects;
    double totalAngleDefect;
    std::vector<double> scalarMeanCurvatures;
    std::vector<double> circumcentricDualAreas;
    std::vector<double> minPrincipalCurvatures;
    std::vector<double> maxPrincipalCurvatures;

    // Constructor
    DiscreteCurvaturesAndNormalsTest() {

        std::vector<Vector3> v;
        std::vector<std::array<int, 3>> f;

        std::string filepath = "../include/test-curv-soln.txt";
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

                } else if (X == "angle") {
                    double val;
                    iss >> val;
                    angles.push_back(val);

                } else if (X == "dihedralAngle") {
                    double val;
                    iss >> val;
                    dihedralAngles.push_back(val);

                } else if (X == "angleWeightedNormal") {
                    double x, y, z;
                    iss >> x >> y >> z;
                    angleWeightedNormals.push_back({x, y, z});

                } else if (X == "sphereInscribedNormal") {
                    double x, y, z;
                    iss >> x >> y >> z;
                    sphereInscribedNormals.push_back({x, y, z});

                } else if (X == "areaWeightedNormal") {
                    double x, y, z;
                    iss >> x >> y >> z;
                    areaWeightedNormals.push_back({x, y, z});

                } else if (X == "gaussCurvatureNormal") {
                    double x, y, z;
                    iss >> x >> y >> z;
                    gaussCurvatureNormals.push_back({x, y, z});

                } else if (X == "meanCurvatureNormal") {
                    double x, y, z;
                    iss >> x >> y >> z;
                    meanCurvatureNormals.push_back({x, y, z});

                } else if (X == "angleDefect") {
                    double val;
                    iss >> val;
                    angleDefects.push_back(val);

                } else if (X == "totalAngleDefect") {
                    double val;
                    iss >> val;
                    totalAngleDefect = val;

                } else if (X == "scalarMeanCurvature") {
                    double val;
                    iss >> val;
                    scalarMeanCurvatures.push_back(val);

                } else if (X == "circumcentricDualArea") {
                    double val;
                    iss >> val;
                    circumcentricDualAreas.push_back(val);

                } else if (X == "k1") {
                    double val;
                    iss >> val;
                    minPrincipalCurvatures.push_back(val);

                } else if (X == "k2") {
                    double val;
                    iss >> val;
                    maxPrincipalCurvatures.push_back(val);
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
        std::tie(this->mesh, this->geometry) = makeManifoldSurfaceMeshAndGeometry(vMat, fMat);
    }

    virtual ~DiscreteCurvaturesAndNormalsTest() {}
};

TEST_F(DiscreteCurvaturesAndNormalsTest, angle) {

    bool allAnglesCorrect = true;
    for (Corner c : mesh->corners()) {
        if (abs(geometry->angle(c) - angles[c.getIndex()]) > 1e-5) {
            allAnglesCorrect = false;
            break;
        }
    }
    EXPECT_TRUE(allAnglesCorrect) << "Angle at a corner is incorrect";
}

TEST_F(DiscreteCurvaturesAndNormalsTest, dihedralAngle) {

    bool allDihedralAnglesCorrect = true;
    for (Halfedge he : mesh->halfedges()) {
        if (abs(geometry->dihedralAngle(he) - dihedralAngles[he.getIndex()]) > 1e-4) {
            allDihedralAnglesCorrect = false;
            break;
        }
    }
    EXPECT_TRUE(allDihedralAnglesCorrect) << "Dihedral angle at a halfedge is incorrect";
}

TEST_F(DiscreteCurvaturesAndNormalsTest, vertexNormalAngleWeighted) {

    bool allNormalsCorrect = true;
    for (Vertex v : mesh->vertices()) {
        if ((geometry->vertexNormalAngleWeighted(v) - angleWeightedNormals[v.getIndex()]).norm() > 1e-5) {
            allNormalsCorrect = false;
            break;
        }
    }
    EXPECT_TRUE(allNormalsCorrect) << "Normal at a vertex computed using the 'tip angle weights' method is incorrect";
}

TEST_F(DiscreteCurvaturesAndNormalsTest, vertexNormalSphereInscribed) {

    bool allNormalsCorrect = true;
    for (Vertex v : mesh->vertices()) {
        if ((geometry->vertexNormalSphereInscribed(v) - sphereInscribedNormals[v.getIndex()]).norm() > 1e-5) {
            allNormalsCorrect = false;
            break;
        }
    }
    EXPECT_TRUE(allNormalsCorrect) << "Normal at a vertex computed using the 'inscribed sphere' method is incorrect";
}

TEST_F(DiscreteCurvaturesAndNormalsTest, vertexNormalAreaWeighted) {

    bool allNormalsCorrect = true;
    for (Vertex v : mesh->vertices()) {
        if ((geometry->vertexNormalAreaWeighted(v) - areaWeightedNormals[v.getIndex()]).norm() > 1e-5) {
            allNormalsCorrect = false;
            break;
        }
    }
    EXPECT_TRUE(allNormalsCorrect) << "Normal at a vertex computed using the 'face area weights' method is incorrect";
}

TEST_F(DiscreteCurvaturesAndNormalsTest, vertexNormalGaussianCurvature) {

    bool allNormalsCorrect = true;
    for (Vertex v : mesh->vertices()) {
        // Allow defining outward or inward-pointing normals
        if ((geometry->vertexNormalGaussianCurvature(v) - gaussCurvatureNormals[v.getIndex()]).norm() > 1e-5 &&
            (geometry->vertexNormalGaussianCurvature(v) + gaussCurvatureNormals[v.getIndex()]).norm() > 1e-5) {
            allNormalsCorrect = false;
            break;
        }
    }
    EXPECT_TRUE(allNormalsCorrect) << "Normal at a vertex computed using the 'Gauss curvature' method is incorrect";
}

TEST_F(DiscreteCurvaturesAndNormalsTest, vertexNormalMeanCurvature) {

    bool allNormalsCorrect = true;
    for (Vertex v : mesh->vertices()) {
        // Allow defining outward or inward-pointing normals
        if ((geometry->vertexNormalMeanCurvature(v) - meanCurvatureNormals[v.getIndex()]).norm() > 1e-5 &&
            (geometry->vertexNormalMeanCurvature(v) + meanCurvatureNormals[v.getIndex()]).norm() > 1e-5) {
            allNormalsCorrect = false;
            break;
        }
    }
    EXPECT_TRUE(allNormalsCorrect) << "Normal at a vertex computed using the 'Mean curvature' method is incorrect";
}

TEST_F(DiscreteCurvaturesAndNormalsTest, angleDefect) {

    bool allAngleDefectsCorrect = true;
    for (Vertex v : mesh->vertices()) {
        if (abs(geometry->angleDefect(v) - angleDefects[v.getIndex()]) > 1e-5) {
            allAngleDefectsCorrect = false;
            break;
        }
    }
    EXPECT_TRUE(allAngleDefectsCorrect) << "Angle defect at a vertex is incorrect";
}

TEST_F(DiscreteCurvaturesAndNormalsTest, totalAngleDefect) {

    EXPECT_TRUE(abs(totalAngleDefect - geometry->totalAngleDefect()) < 1e-5) << "Total angle defect is incorrect";
}

TEST_F(DiscreteCurvaturesAndNormalsTest, scalarMeanCurvature) {

    bool allCurvaturesCorrect = true;
    for (Vertex v : mesh->vertices()) {
        if (abs(geometry->scalarMeanCurvature(v) - scalarMeanCurvatures[v.getIndex()]) > 1e-5) {
            allCurvaturesCorrect = false;
            break;
        }
    }
    EXPECT_TRUE(allCurvaturesCorrect) << "Scalar mean curvature at a vertex is incorrect";
}

TEST_F(DiscreteCurvaturesAndNormalsTest, circumcentricDualArea) {

    bool allCircumcentricDualAreasCorrect = true;
    for (Vertex v : mesh->vertices()) {
        if (abs(geometry->circumcentricDualArea(v) - circumcentricDualAreas[v.getIndex()]) > 1e-5) {
            allCircumcentricDualAreasCorrect = false;
            break;
        }
    }
    EXPECT_TRUE(allCircumcentricDualAreasCorrect) << "Circumcentric dual area of a vertex is incorrect";
}

TEST_F(DiscreteCurvaturesAndNormalsTest, principalCurvatures) {

    bool allKMinCorrect = true;
    bool allKMaxCorrect = true;
    for (Vertex v : mesh->vertices()) {
        std::pair<double, double> K = geometry->principalCurvatures(v);
        double k1 = std::min(K.first, K.second);
        double k2 = std::max(K.first, K.second);
        double a = minPrincipalCurvatures[v.getIndex()];
        double b = maxPrincipalCurvatures[v.getIndex()];
        double k1_soln = std::min(a, b);
        double k2_soln = std::max(a, b);
        if (abs(k1 - k1_soln) > 1e-5) {
            allKMinCorrect = false;
        }
        if (abs(k2 - k2_soln) > 1e-5) {
            allKMaxCorrect = false;
        }
        if (!allKMinCorrect && !allKMaxCorrect) {
            break;
        }
    }
    EXPECT_TRUE(allKMinCorrect) << "Min principal curvature at a vertex is incorrect";
    EXPECT_TRUE(allKMaxCorrect) << "Max principal curvature at a vertex is incorrect";
}
} // namespace


int main(int argc, char** argv) {

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}