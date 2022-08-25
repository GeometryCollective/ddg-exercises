// POISSON PROBLEM

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "colormap.h"
#include "scalar-poisson-problem.h"
#include <algorithm>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;
// so we can more easily pass these to different classes
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;
std::string MESHNAME;

// Some global variables
ScalarPoissonProblem SPP;
Vector<double> RHO;
polyscope::SurfaceVertexColorQuantity* solnColors;
double DENSITY = 1.0;                      // essentially a dummy variable
polyscope::SurfaceGraphQuantity* currVert; // currently active vertex
double vertexRadius;


glm::vec<3, float> ORANGE = {1, 0.65, 0};


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
 * Show selected vertices.
 * This function gets called every time an element is selected on-screen.
 */
void showSelected() {

    // If we've de-selected a vertex, reset the corresponding entry in <rho>.
    int toDelete = polyscope::state::deleteVertexIndex;
    if (toDelete != -1) {
        RHO[(size_t)toDelete] = 0;
    }

    // Show selected vertices in yellow
    std::vector<Vector3> vertPos;
    std::vector<std::array<size_t, 2>> vertInd;
    for (std::set<size_t>::iterator it = polyscope::state::subset.vertices.begin();
         it != polyscope::state::subset.vertices.end(); ++it) {
        vertPos.push_back(geometry->inputVertexPositions[*it]);
        // default density is 1
        if (RHO[*it] == 0) {
            RHO[*it] = 1;
        }
    }
    polyscope::SurfaceGraphQuantity* showVerts = psMesh->addSurfaceGraphQuantity("selected vertices", vertPos, vertInd);
    showVerts->setEnabled(true);
    showVerts->setRadius(vertexRadius);
    showVerts->setColor(ORANGE);

    // Show the currently selected vertex in red.
    int currIdx = polyscope::state::currVertexIndex;
    if (currIdx != -1) {
        std::vector<Vector3> pos = {geometry->inputVertexPositions[currIdx]};
        currVert = psMesh->addSurfaceGraphQuantity("current vertex", pos, std::vector<std::array<size_t, 2>>());
        currVert->setEnabled(true);
        currVert->setRadius(vertexRadius);
        currVert->setColor({1.0, 0.0, 0.0});
    } else {
        currVert->setEnabled(false);
    }
}

void redraw() {
    showSelected();
    polyscope::requestRedraw();
}

/*
 * Map solution to mesh colors.
 */
std::vector<std::array<double, 3>> computeColors(const Vector<double>& sol) {

    // Determine maximum-magnitude element for scaling purposes
    double maxVal = 0;
    for (size_t i = 0; i < mesh->nVertices(); i++) {
        maxVal = std::max(maxVal, abs(sol[i]));
    }
    if (maxVal < 1e-3) {
        maxVal = 1e-3;
    }

    std::vector<std::array<double, 3>> colors;
    for (Vertex v : mesh->vertices()) {
        colors.push_back(mapToColor(sol[v.getIndex()], -maxVal, maxVal, "seismic"));
    }
    return colors;
}

// Set mesh color.
void setColors(const std::vector<std::array<double, 3>>& colors) {
    solnColors = psMesh->addVertexColorQuantity("Solution", colors);
    solnColors->setEnabled(true);
}

// Solve and update mesh colors.
void update() {
    if (polyscope::state::subset.vertices.size() > 0) {
        setColors(computeColors(SPP.solve(RHO)));
    }
}

void functionCallback() {

    int idx = polyscope::state::currVertexIndex;
    ImGui::Text("Density: ");
    ImGui::SameLine(); // put label on the left
    ImGui::InputDouble("", idx == -1 ? &DENSITY : &RHO[idx], -std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity(), "%0.2f");

    // Solve and update colors
    if (ImGui::Button("Solve")) {
        update();
        redraw();
    }
    if (ImGui::Button("Reset")) {
        polyscope::state::subset.vertices.clear();
        polyscope::state::currVertexIndex = -1;
        RHO = Vector<double>::Zero(mesh->nVertices());
        psMesh->setSurfaceColor({1.0, 1.0, 1.0});
        solnColors->setEnabled(false);
        redraw();
    }
}


int main(int argc, char** argv) {
    // Configure the argument parser
    args::ArgumentParser parser("15-458 HW3");
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
    std::string filepath = "../../../input/bunny.obj";
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

    // Add visualization options.
    flipZ();
    vertexRadius = geometry->meanEdgeLength() * 0.2;
    psMesh->setSmoothShade(true);
    psMesh->setSurfaceColor({1.0, 1.0, 1.0});
    // initialize to something
    currVert =
        psMesh->addSurfaceGraphQuantity("current vertex", std::vector<Vector3>(), std::vector<std::array<size_t, 2>>());
    // initialize to something in case "Reset" is pressed before anything happens
    solnColors = psMesh->addVertexColorQuantity("Solution", std::vector<std::array<double, 3>>(mesh->nVertices()));
    SPP = ScalarPoissonProblem(mesh, geometry);
    RHO = Vector<double>::Zero(mesh->nVertices());

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}