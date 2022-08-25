// HEAT METHOD

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "colormap.h"
#include "heat-method.h"

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
HeatMethod HM;
Vector<double> DELTA;                      // sources
polyscope::SurfaceGraphQuantity* currVert; // currently active vertex
Vector<double> SOLUTION;
polyscope::SurfaceVertexColorQuantity* solnColors;
polyscope::SurfaceGraphQuantity* isolines;
double maxPhi = 0.0;
double vertexRadius;
double isolinesRadius;


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
 * Map solution to mesh colors.
 */
std::vector<std::array<double, 3>> computeColors(const Vector<double>& sol) {

    // Determine maximum-magnitude element for scaling purposes
    maxPhi = 0;
    for (size_t i = 0; i < mesh->nVertices(); i++) {
        maxPhi = std::max(maxPhi, sol[i]); // distances should be >= 0
    }

    std::vector<std::array<double, 3>> colors;
    for (Vertex v : mesh->vertices()) {
        colors.push_back(mapToColor(maxPhi - sol[v.getIndex()], 0, maxPhi, "hot")); // invert colormap
    }
    return colors;
}

// Set mesh color.
void setColors(const std::vector<std::array<double, 3>>& colors) {
    solnColors = psMesh->addVertexColorQuantity("Solution", colors);
    solnColors->setEnabled(true);
}


/*
 * Display isolines.
 */
void showIsolines() {

    std::vector<Vector3> positions;
    std::vector<std::array<size_t, 2>> edgeInds;
    double distBetweenLines = maxPhi / 20.0; // enforce spacing
    for (Face f : mesh->faces()) {
        std::vector<Vector3> pos;
        for (Halfedge he : f.adjacentHalfedges()) {
            double vs = SOLUTION[he.tailVertex().getIndex()];
            double vd = SOLUTION[he.tipVertex().getIndex()];
            int region1 = floor(vs / distBetweenLines);
            int region2 = floor(vd / distBetweenLines);
            if (region1 != region2) {
                double val = region2 * distBetweenLines;
                if (region1 > region2) {
                    val = region1 * distBetweenLines;
                }
                double t = (val - vs) / (vd - vs);
                Vector3 ps = geometry->inputVertexPositions[he.tailVertex()];
                Vector3 pd = geometry->inputVertexPositions[he.tipVertex()];
                Vector3 p = ps + t * (pd - ps);
                pos.push_back(p);
            }
        }
        if (pos.size() == 2) {
            positions.push_back(pos[0]);
            positions.push_back(pos[1]);
            edgeInds.push_back({positions.size() - 2, positions.size() - 1});
        }
    }
    isolines = psMesh->addSurfaceGraphQuantity("Isolines", positions, edgeInds);
    isolines->setEnabled(true);
    isolines->setRadius(isolinesRadius);
    isolines->setColor({0.0, 0.0, 0.0});
}

/*
 * Show selected vertices.
 * This function gets called every time an element is selected on-screen.
 */
void showSelected() {

    // Show selected vertices in yellow
    std::vector<Vector3> vertPos;
    std::vector<std::array<size_t, 2>> vertInd;
    DELTA = Vector<double>::Zero(mesh->nVertices());
    for (std::set<size_t>::iterator it = polyscope::state::subset.vertices.begin();
         it != polyscope::state::subset.vertices.end(); ++it) {
        vertPos.push_back(geometry->inputVertexPositions[*it]);
        DELTA[*it] = 1;
    }
    polyscope::SurfaceGraphQuantity* showVerts = psMesh->addSurfaceGraphQuantity("selected vertices", vertPos, vertInd);
    showVerts->setEnabled(true);
    showVerts->setRadius(vertexRadius);
    showVerts->setColor({1.0, 0.65, 0.0});

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

void functionCallback() {

    if (ImGui::Button("Solve")) {
        if (polyscope::state::subset.vertices.size() > 0) {
            SOLUTION = HM.compute(DELTA);
            if (SOLUTION.norm() > 0) {
                setColors(computeColors(SOLUTION));
                showIsolines();
                redraw();
            }
        }
    }
    if (ImGui::Button("Reset")) {
        polyscope::state::subset.vertices.clear();
        polyscope::state::currVertexIndex = -1;
        psMesh->setSurfaceColor({1.0, 0.45, 0.0});
        solnColors->setEnabled(false);
        isolines->setEnabled(false);
        redraw();
    }
}


int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("15-458 HW5");
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

    // Initalize
    flipZ();
    double lengthScale = geometry->meanEdgeLength();
    vertexRadius = lengthScale * 0.2;
    isolinesRadius = lengthScale * 0.05;
    psMesh->setSmoothShade(true);
    psMesh->setSurfaceColor({1.0, 0.45, 0.0}); // orange
    currVert =
        psMesh->addSurfaceGraphQuantity("current vertex", std::vector<Vector3>(), std::vector<std::array<size_t, 2>>());
    // initialize to something in case "Reset" is pressed before anything happens
    solnColors = psMesh->addVertexColorQuantity("Solution", std::vector<std::array<double, 3>>(mesh->nVertices()));
    isolines =
        psMesh->addSurfaceGraphQuantity("Isolines", std::vector<Vector3>(), std::vector<std::array<size_t, 2>>());
    HM = HeatMethod(mesh, geometry);
    DELTA = Vector<double>::Zero(mesh->nVertices());

    // Give control to the polyscope gui
    polyscope::show();

    delete mesh;
    delete geometry;

    return EXIT_SUCCESS;
}