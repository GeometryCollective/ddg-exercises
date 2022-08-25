// SIMPLICIAL COMPLEX

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;
// so we can more easily pass these to different classes while preserving syntax
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;
std::string MESHNAME;

// Some global variables
SimplicialComplexOperators SCO;
bool isComplexResult = false;
int isPureComplexResult = -1;
double vertexRadius;
double edgeRadius;


std::array<double, 3> BLUE = {0.11, 0.388, 0.89};
glm::vec<3, float> ORANGE_VEC = {1, 0.65, 0};
std::array<double, 3> ORANGE = {1, 0.65, 0};


void flipZ() {
    // Rotate mesh 180 deg about up-axis on startup
    glm::mat4x4 rot = glm::rotate(glm::mat4x4(1.0f), static_cast<float>(PI), glm::vec3(0, 1, 0));
    for (Vertex v : mesh->vertices()) {
        Vector3 vec = geometry->inputVertexPositions[v];
        glm::vec4 rvec = {vec[0], vec[1], vec[2], 1.0};
        rvec = rot * rvec;
        geometry->inputVertexPositions[v] = {rvec[0], rvec[1], rvec[2]};
    }
    psMesh->updateVertexPositions(geometry->inputVertexPositions);
}

/*
 * Display the selected simplices.
 * TODO: Use SurfaceVertexCountQuantity* SurfaceMesh::addVertexCountQuantity, etc. instead of SurfaceGraphQuantity for
 * cleaner code
 */
void showSelected() {

    // Show selected vertices.
    std::vector<Vector3> vertPos;
    std::vector<std::array<size_t, 2>> vertInd;
    for (std::set<size_t>::iterator it = polyscope::state::subset.vertices.begin();
         it != polyscope::state::subset.vertices.end(); ++it) {
        vertPos.push_back(geometry->inputVertexPositions[*it]);
    }
    polyscope::SurfaceGraphQuantity* showVerts = psMesh->addSurfaceGraphQuantity("selected vertices", vertPos, vertInd);
    showVerts->setEnabled(true);
    showVerts->setRadius(vertexRadius);
    showVerts->setColor(ORANGE_VEC);

    // Show selected edges.
    std::vector<Vector3> edgePos;
    std::vector<std::array<size_t, 2>> edgeInd;
    for (std::set<size_t>::iterator it = polyscope::state::subset.edges.begin();
         it != polyscope::state::subset.edges.end(); ++it) {
        Edge e = mesh->edge(*it);
        edgePos.push_back(geometry->inputVertexPositions[e.firstVertex()]);
        edgePos.push_back(geometry->inputVertexPositions[e.secondVertex()]);
        size_t i = edgeInd.size();
        edgeInd.push_back({2 * i, 2 * i + 1});
    }
    polyscope::SurfaceGraphQuantity* showEdges = psMesh->addSurfaceGraphQuantity("selected edges", edgePos, edgeInd);
    showEdges->setEnabled(true);
    showEdges->setRadius(edgeRadius);
    showEdges->setColor(ORANGE_VEC);

    // Show selected faces.
    std::vector<std::array<double, 3>> faceColors(mesh->nFaces());
    for (size_t i = 0; i < mesh->nFaces(); i++) {
        faceColors[i] = BLUE;
    }
    for (std::set<size_t>::iterator it = polyscope::state::subset.faces.begin();
         it != polyscope::state::subset.faces.end(); ++it) {
        faceColors[*it] = ORANGE;
    }
    polyscope::SurfaceFaceColorQuantity* showFaces = psMesh->addFaceColorQuantity("selected faces", faceColors);
    showFaces->setEnabled(true);
}

void redraw() {
    showSelected();
    polyscope::requestRedraw();
}

/*
 * Buttons for the simplicial operators.
 */
void functionCallback() {

    if (ImGui::Button("Reset")) {
        polyscope::state::subset.vertices.clear();
        polyscope::state::subset.edges.clear();
        polyscope::state::subset.faces.clear();
        redraw();
    }

    if (ImGui::Button("isComplex")) {
        isComplexResult = SCO.isComplex(polyscope::state::subset);
    }
    ImGui::SameLine(100);
    ImGui::Text(isComplexResult ? "true" : "false");

    if (ImGui::Button("isPureComplex")) {
        isPureComplexResult = SCO.isPureComplex(polyscope::state::subset);
    }
    ImGui::SameLine(130);
    ImGui::Text("%d", isPureComplexResult);

    if (ImGui::Button("Boundary")) {
        MeshSubset S = SCO.boundary(polyscope::state::subset);
        polyscope::state::subset = S;
        redraw();
    }

    if (ImGui::Button("Star")) {
        MeshSubset S = SCO.star(polyscope::state::subset);
        polyscope::state::subset = S;
        redraw();
    }

    if (ImGui::Button("Closure")) {
        MeshSubset S = SCO.closure(polyscope::state::subset);
        polyscope::state::subset = S;
        redraw();
    }

    if (ImGui::Button("Link")) {
        MeshSubset S = SCO.link(polyscope::state::subset);
        polyscope::state::subset = S;
        redraw();
    }
}


int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("15-458 HW0");
    args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help&) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    // If a mesh name was not given, use default mesh.
    std::string filepath = "../../../input/small_disk.obj";
    if (inputFilename) {
        filepath = args::get(inputFilename);
    }

    MESHNAME = polyscope::guessNiceNameFromPath(filepath);

    // Load mesh
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();

    // Get indices for element picking
    polyscope::state::facePickIndStart = mesh->nVertices();
    polyscope::state::edgePickIndStart = polyscope::state::facePickIndStart + mesh->nFaces();
    polyscope::state::halfedgePickIndStart = polyscope::state::edgePickIndStart + mesh->nEdges();

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = functionCallback;

    // Add mesh to GUI
    psMesh = polyscope::registerSurfaceMesh(MESHNAME, geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));

    // Mesh initialization
    SCO.initialize(mesh, geometry);

    // Add visualization options.
    flipZ();
    double lengthScale = geometry->meanEdgeLength();
    vertexRadius = lengthScale * 0.1;
    edgeRadius = lengthScale * 0.05;

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}