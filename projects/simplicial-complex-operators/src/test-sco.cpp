#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "simplicial-complex-operators.h"
#include "gtest/gtest.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace {

class SimplicialComplexOperatorsTest : public ::testing::Test {

  public:
    SimplicialComplexOperators SCO;

  protected:
    /*
     * Load in the bunny mesh to test on.
     */
    SimplicialComplexOperatorsTest() {

        std::string filepath = "../../../input/bunny.obj";
        std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
        std::unique_ptr<VertexPositionGeometry> geometry_uptr;
        std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
        ManifoldSurfaceMesh* mesh = mesh_uptr.release();
        VertexPositionGeometry* geometry = geometry_uptr.release();
        SCO.initialize(mesh, geometry);

        // In case students are accessing using cached arrays.
        SCO.geometry->requireVertexIndices();
        SCO.geometry->requireEdgeIndices();
        SCO.geometry->requireFaceIndices();
    }

    virtual ~SimplicialComplexOperatorsTest() {}

    /*
     * Construct a MeshSubset with the given vertices, edges, and faces, to test the given function on. Check if we get
     * the expected result. Designed to test isComplex() and isPureComplex().
     */
    void testOnMeshSubset(const std::string& functionName, const std::set<size_t>& vertices,
                          const std::set<size_t>& edges, const std::set<size_t>& faces, int expectedResult = -1,
                          const std::string subcomplexDescription = "") {

        // Student indexing should be same as one automically assigned by Geometry Central; students don't have to do
        // that function in the C++ version.
        MeshSubset S = MeshSubset(vertices, edges, faces);

        std::cerr << "\t A" << subcomplexDescription << ":\n";
        if (functionName == "isComplex") {
            EXPECT_EQ(expectedResult, SCO.isComplex(S))
                << "\t\t isComplex() returns wrong result for a" << subcomplexDescription;
        } else if (functionName == "isPureComplex") {
            EXPECT_EQ(expectedResult, SCO.isPureComplex(S))
                << "\t\t isPureComplex() returns wrong result for a" << subcomplexDescription;
        }
    }

    /*
     * Pretty much the same as testOnMeshSubset(), but specifically customized for testing star(), closure(), link(),
     * and boundary(), which return MeshSubsets.
     */
    void testReturnedMeshSubset(const std::string& functionName, const std::set<size_t>& vertices,
                                const std::set<size_t>& edges, const std::set<size_t>& faces,
                                const std::set<size_t>& expectedVertices, const std::set<size_t>& expectedEdges,
                                const std::set<size_t>& expectedFaces) {

        // Student indexing should be same as one automically assigned by Geometry Central; students don't have to do
        // that function in the C++ version.
        MeshSubset S = MeshSubset(vertices, edges, faces);
        MeshSubset result;
        if (functionName == "star") {
            result = SCO.star(S);
        } else if (functionName == "link") {
            result = SCO.link(S);
        } else if (functionName == "closure") {
            result = SCO.closure(S);
        } else if (functionName == "boundary") {
            result = SCO.boundary(S);
        }

        bool verticesAreCorrect = result.vertices == expectedVertices;
        bool edgesAreCorrect = result.edges == expectedEdges;
        bool facesAreCorrect = result.faces == expectedFaces;

        EXPECT_TRUE(verticesAreCorrect) << "\t\t The vertices in your " << functionName << " are wrong";
        EXPECT_TRUE(edgesAreCorrect) << "\t\t The edges in your " << functionName << " are wrong";
        EXPECT_TRUE(facesAreCorrect) << "\t\t The faces in your " << functionName << " are wrong";
    }
};

TEST_F(SimplicialComplexOperatorsTest, isComplex) {

    std::cerr << "Testing isComplex()..." << std::endl;

    testOnMeshSubset("isComplex", {615}, {}, {}, true, " vertex");
    testOnMeshSubset("isComplex", {}, {607}, {}, false, "n edge");
    testOnMeshSubset("isComplex", {8880, 8996}, {1988}, {}, true, " closed edge");
    testOnMeshSubset("isComplex", {}, {}, {584}, false, " face");
    testOnMeshSubset("isComplex", {}, {1988, 1989, 1987}, {711}, false, " face with its edges");
    testOnMeshSubset("isComplex", {8880, 8620, 8996}, {1988, 1989, 1987}, {711}, true, " closed face");
}

TEST_F(SimplicialComplexOperatorsTest, isPureComplex) {

    std::cerr << "Testing isPureComplex()..." << std::endl;

    testOnMeshSubset("isPureComplex", {731}, {}, {}, 0, " vertex");
    testOnMeshSubset("isPureComplex", {}, {1969}, {}, -1, "n edge");
    testOnMeshSubset("isPureComplex", {9298, 2450}, {6647}, {}, 1, " closed edge");
    testOnMeshSubset("isPureComplex", {}, {}, {1261}, -1, " face");
    testOnMeshSubset("isPureComplex", {}, {1988, 1989, 1987}, {711}, -1, " face with its edges");
    testOnMeshSubset("isPureComplex", {8880, 8620, 8996}, {1988, 1989, 1987}, {711}, 2, " closed face");
    testOnMeshSubset("isPureComplex", {8880, 8620, 8996, 9298, 2450}, {1988, 1989, 1987, 6647}, {711}, -1,
                     " closed face and closed edge");
}

TEST_F(SimplicialComplexOperatorsTest, A0) {

    std::cerr << "Testing buildVertexEdgeAdjacencyMatrix()..." << std::endl;

    std::cerr << "\t Has |E| rows:" << std::endl;
    EXPECT_EQ(SCO.A0.rows(), SCO.mesh->nEdges()) << "\t\t A0 does not have |E| rows";

    std::cerr << "\t Has |V| columns:" << std::endl;
    EXPECT_EQ(SCO.A0.cols(), SCO.mesh->nVertices()) << "\t\t A0 does not have |V| columns";

    std::cerr << "\t Rows sum to 2:" << std::endl;
    Vector<size_t> rowSums = SCO.A0 * Vector<size_t>::Ones(SCO.mesh->nVertices());
    bool rowSumsAllTwo = true;
    for (int i = 0; i < SCO.A0.rows(); i++) {
        if (rowSums[i] != 2) {
            rowSumsAllTwo = false;
            break;
        }
    }
    EXPECT_TRUE(rowSumsAllTwo) << "\t\t There is a row in A0 which does not sum to 2.";

    std::cerr << "\t Columns sum to degrees:" << std::endl;
    Vector<size_t> colSums = (Vector<size_t>::Ones(SCO.mesh->nEdges()).transpose() * SCO.A0);
    bool colSumsAllCorrect = true;
    for (int i = 0; i < SCO.A0.cols(); i++) {
        if (colSums[i] != SCO.mesh->vertex(i).degree()) {
            colSumsAllCorrect = false;
            break;
        }
    }
    EXPECT_TRUE(colSumsAllCorrect)
        << "\t\t There is a column in your A0 which does not sum to the degree of the corresponding vertex.";
}

TEST_F(SimplicialComplexOperatorsTest, A1) {

    std::cerr << "Testing buildFaceEdgeAdjacencyMatrix()..." << std::endl;

    std::cerr << "\t Has |F| rows:" << std::endl;
    EXPECT_EQ(SCO.A1.rows(), SCO.mesh->nFaces()) << "\t\t A1 does not have |F| rows";

    std::cerr << "\t Has |E| columns:" << std::endl;
    EXPECT_EQ(SCO.A1.cols(), SCO.mesh->nEdges()) << "\t\t A1 does not have |E| columns";

    std::cerr << "\t Rows sum to number of faces:" << std::endl;
    Vector<size_t> rowSums = SCO.A1 * Vector<size_t>::Ones(SCO.mesh->nEdges());
    bool rowSumsAllCorrect = true;
    for (int i = 0; i < SCO.A1.rows(); i++) {
        if (rowSums[i] != SCO.mesh->face(i).degree()) {
            rowSumsAllCorrect = false;
            break;
        }
    }
    EXPECT_TRUE(rowSumsAllCorrect)
        << "\t\t There is a row in A1 which does not sum to the number of edges of the corresponding face";

    std::cerr << "\t Columns sum to 2:" << std::endl;
    Vector<size_t> colSums = (Vector<size_t>::Ones(SCO.mesh->nFaces()).transpose() * SCO.A1);
    bool colSumsAllCorrect = true;
    for (int i = 0; i < SCO.A1.cols(); i++) {
        if (colSums[i] != 2) {
            colSumsAllCorrect = false;
            break;
        }
    }
    EXPECT_TRUE(colSumsAllCorrect) << "\t\t There is a column in A1 which does not sum to 2";
}

TEST_F(SimplicialComplexOperatorsTest, buildVertexVector) {

    std::cerr << "Testing buildVertexVector()..." << std::endl;
    MeshSubset S = MeshSubset({1, 2, 3}, {1, 2, 3}, {1, 2, 3});

    // Could add more cases to testOnMeshSubset(), but after a certain point it gets inefficient
    std::cerr << "\t Correct size:" << std::endl;
    Vector<size_t> vec = SCO.buildVertexVector(S);
    EXPECT_EQ(SCO.mesh->nVertices(), vec.rows()) << "\t\t Vertex vector has wrong size";

    std::cerr << "\t All entries 0 or 1:" << std::endl;
    bool entriesAreBinary = true;
    for (int i = 0; i < vec.rows(); i++) {
        if (vec[i] != 1 && vec[i] != 0) {
            entriesAreBinary = false;
            break;
        }
    }
    EXPECT_TRUE(entriesAreBinary) << "\t\t Vertex vector has entry that is not 0 or 1";
}

TEST_F(SimplicialComplexOperatorsTest, buildEdgeVector) {

    std::cerr << "Testing buildEdgeVector()..." << std::endl;
    MeshSubset S = MeshSubset({1, 2, 3}, {1, 2, 3}, {1, 2, 3});

    std::cerr << "\t Correct size:" << std::endl;
    Vector<size_t> vec = SCO.buildEdgeVector(S);
    EXPECT_EQ(SCO.mesh->nEdges(), vec.rows()) << "\t\t Edge vector has wrong size";

    std::cerr << "\t All entries 0 or 1:" << std::endl;
    bool entriesAreBinary = true;
    for (int i = 0; i < vec.rows(); i++) {
        if (vec[i] != 1 && vec[i] != 0) {
            entriesAreBinary = false;
            break;
        }
    }
    EXPECT_TRUE(entriesAreBinary) << "\t\t Edge vector has entry that is not 0 or 1";
}

TEST_F(SimplicialComplexOperatorsTest, buildFaceVector) {

    std::cerr << "Testing buildFaceVector()..." << std::endl;
    MeshSubset S = MeshSubset({1, 2, 3}, {1, 2, 3}, {1, 2, 3});

    std::cerr << "\t Correct size:" << std::endl;
    Vector<size_t> vec = SCO.buildFaceVector(S);
    EXPECT_EQ(SCO.mesh->nFaces(), vec.rows()) << "\t\t Face vector has wrong size";

    std::cerr << "\t All entries 0 or 1:" << std::endl;
    bool entriesAreBinary = true;
    for (int i = 0; i < vec.rows(); i++) {
        if (vec[i] != 1 && vec[i] != 0) {
            entriesAreBinary = false;
            break;
        }
    }
    EXPECT_TRUE(entriesAreBinary) << "\t\t Face vector has entry that is not 0 or 1";
}

TEST_F(SimplicialComplexOperatorsTest, boundary) {

    std::cerr << "Testing boundary()..." << std::endl;
    std::set<size_t> vertices = {2026, 3598, 3791, 3792, 8480, 8615, 8729, 8737, 8748,  8866,
                                 8880, 8894, 8996, 9015, 9113, 9172, 9185, 9313, 10054, 12147};
    std::set<size_t> edges = {1988,  4246,  4247,  4248,  6741,  7138,  7139,  7140,  7865,  13648, 15454,
                              15455, 15456, 21304, 21305, 21306, 21376, 21377, 21378, 21587, 21588, 21589,
                              21590, 21591, 21613, 21614, 21615, 22168, 22189, 22190, 22473, 22520, 22529,
                              22530, 23358, 23501, 23502, 24999, 26859, 26860, 26861, 27971, 28011, 29964};
    std::set<size_t> faces = {1597,  2797,  6604,  9512,  9550,  9663,  9664,  9675,  9820,  9872,  9993,  10147, 10175,
                              10619, 10693, 11524, 11617, 12178, 12716, 13406, 13438, 14622, 16029, 17736, 17884};
    MeshSubset S = MeshSubset(vertices, edges, faces);

    std::cerr << "\t Squares to zero:" << std::endl;
    MeshSubset result = SCO.boundary(SCO.boundary(S));
    EXPECT_EQ(result.vertices.size(), 0) << "Double boundary of a set contains extra vertices";
    EXPECT_EQ(result.edges.size(), 0) << "Double boundary of a set contains extra edges";
    EXPECT_EQ(result.faces.size(), 0) << "Double boundary of a set contains extra faces";

    std::cerr << "\t Boundary of a vertex:" << std::endl;
    testReturnedMeshSubset("boundary", {915}, {}, {}, {}, {}, {});
    std::cerr << "\t Boundary of an edge:" << std::endl;
    testReturnedMeshSubset("boundary", {8621, 8620}, {6381}, {}, {8621, 8620}, {}, {});
    std::cerr << "\t Boundary of a face:" << std::endl;
    testReturnedMeshSubset("boundary", {9484, 9483, 9112}, {9384, 9385, 9386}, {3744}, {9484, 9483, 9112},
                           {9384, 9385, 9386}, {});
    std::cerr << "\t Boundary of adjacent faces:" << std::endl;
    testReturnedMeshSubset("boundary", {8669, 12283, 12306, 12308, 12347},
                           {14904, 15159, 22231, 23221, 32126, 32252, 35098}, {15938, 18115, 25491},
                           {8669, 12283, 12306, 12308, 12347}, {14904, 15159, 22231, 23221, 32252}, {});
}

TEST_F(SimplicialComplexOperatorsTest, star) {

    std::cerr << "Testing star()..." << std::endl;
    std::cerr << "\t Star of a vertex:" << std::endl;
    testReturnedMeshSubset("star", {12311}, {}, {}, {12311}, {15481, 15482, 21057, 21529, 21530, 21586},
                           {6619, 9381, 9638, 9662, 10652, 19575});
    std::cerr << "\t Star of two vertices:" << std::endl;
    testReturnedMeshSubset("star", {12311, 506}, {}, {}, {506, 12311},
                           {15481, 15482, 15516, 15517, 15599, 21057, 21529, 21530, 21531, 21586},
                           {6619, 6635, 6678, 9381, 9638, 9662, 10652, 10742, 19575});
    std::cerr << "\t Star of an edge:" << std::endl;
    testReturnedMeshSubset("star", {}, {21531}, {}, {}, {21531}, {9638, 10742});
    std::cerr << "\t Star of a face:" << std::endl;
    testReturnedMeshSubset("star", {}, {}, {841}, {}, {}, {841});
    std::cerr << "\t Squares to itself:" << std::endl;
    MeshSubset S = MeshSubset({1, 2, 3}, {1, 2, 3}, {1, 2, 3});
    MeshSubset A = SCO.star(S);
    MeshSubset B = SCO.star(A);
    EXPECT_TRUE(A.equals(B)) << "\t\t You star does not square to itself";
}

TEST_F(SimplicialComplexOperatorsTest, closure) {

    std::cerr << "Testing closure()..." << std::endl;
    std::cerr << "\t Closure of a vertex:" << std::endl;
    testReturnedMeshSubset("closure", {794}, {}, {}, {794}, {}, {});
    std::cerr << "\t Closure of two vertices:" << std::endl;
    testReturnedMeshSubset("closure", {731, 1064}, {}, {}, {731, 1064}, {}, {});
    std::cerr << "\t Closure of an edge:" << std::endl;
    testReturnedMeshSubset("closure", {}, {25907}, {}, {12426, 12495}, {25907}, {});
    std::cerr << "\t Closure of a face:" << std::endl;
    testReturnedMeshSubset("closure", {}, {}, {12935}, {12415, 12456, 12540}, {27226, 27227, 27228}, {12935});
    std::cerr << "\t Squares to itself:" << std::endl;
    MeshSubset S = MeshSubset({1, 2, 3}, {1, 2, 3}, {1, 2, 3});
    MeshSubset A = SCO.closure(S);
    MeshSubset B = SCO.closure(A);
    EXPECT_TRUE(A.equals(B)) << "\t\t You closure does not square to itself";
}

TEST_F(SimplicialComplexOperatorsTest, link) {

    std::cerr << "Testing link()..." << std::endl;
    std::cerr << "\t Link of a vertex:" << std::endl;
    testReturnedMeshSubset("link", {12415}, {}, {}, {12347, 12381, 12384, 12456, 12458, 12540},
                           {14903, 21875, 23607, 27226, 27910, 30078}, {});
    std::cerr << "\t Link of two vertices:" << std::endl;
    testReturnedMeshSubset("link", {12426, 12493}, {}, {}, {506, 12391, 12424, 12427, 12492, 12495, 12543, 12565},
                           {3612, 12664, 15517, 15599, 18579, 18698, 21193, 21667}, {});
    std::cerr << "\t Link of an edge:" << std::endl;
    testReturnedMeshSubset("link", {}, {21355}, {}, {12205, 12308}, {}, {});
    std::cerr << "\t Link of a face:" << std::endl;
    testReturnedMeshSubset("link", {}, {}, {1261}, {}, {}, {});
}

} // namespace


int main(int argc, char** argv) {

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}