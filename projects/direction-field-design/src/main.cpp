// DIRECTION FIELD DESIGN

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "trivial-connections.h"

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;
// so we can more easily pass these to different classes
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

// Some global variables
TrivialConnections connections;
Vector<double> RHO;                        // encodes singularity indices
double INDEX = 1.0;                        // for display
polyscope::SurfaceGraphQuantity* currVert; // currently active vertex
double vertexRadius;
polyscope::SurfaceVectorQuantity* directionField;

int EULER_CHARACTERISTIC;
double lengthScale;


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
        // default singularity index is 1
        if (RHO[*it] == 0) {
            RHO[*it] = 1;
        }
    }
    polyscope::SurfaceGraphQuantity* showVerts = psMesh->addSurfaceGraphQuantity("selected vertices", vertPos, vertInd);
    showVerts->setEnabled(true);
    showVerts->setRadius(vertexRadius);
    showVerts->setColor({1, 0.65, 0});

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

void updateDirectionField(const Vector<double>& phi, bool display) {

    std::map<Face, Face> faceParent;
    for (Face f : mesh->faces()) {
        faceParent[f] = f;
    }
    Face root = mesh->face(0);
    FaceData<double> alpha(*mesh);
    alpha[root] = 0;
    std::vector<Face> bag;
    bag.push_back(root);
    while (bag.size() != 0) {
        Face f = bag.back();
        bag.pop_back();
        for (Halfedge he : f.adjacentHalfedges()) {
            Face g = he.twin().face();
            if (faceParent[g] == g && g != root) {
                double sign = he.edge().halfedge() == he ? 1 : -1;
                faceParent[g] = f;
                double connection = sign * phi[he.edge().getIndex()];
                alpha[g] = connections.transportNoRotation(he, alpha[f] - connection);
                bag.push_back(g);
            }
        }
    }

    std::vector<Vector3> field(mesh->nFaces());
    for (Face f : mesh->faces()) {
        Vector3 u = {cos(alpha[f]), sin(alpha[f]), 0};
        Vector3 e1 = geometry->halfedgeVector(f.halfedge()).normalize();
        Vector3 e2 = cross(geometry->faceNormal(f), e1);
        field[f.getIndex()] = e1 * u[0] + e2 * u[1];
    }

    directionField = psMesh->addFaceVectorQuantity("Direction field", field);
    directionField->setVectorColor({0, 0, 0});
    directionField->setVectorLengthScale(lengthScale * 0.5);
    directionField->setEnabled(display);
}

void redraw() {
    showSelected();
    polyscope::requestRedraw();
}

/*
 * User-defined buttons
 */
void functionCallback() {

    ImGui::Text("Euler characteristic: %d", EULER_CHARACTERISTIC);
    ImGui::Text("Sum of singularity indices: %0.2f", RHO.sum());
    ImGui::Text("");

    int idx = polyscope::state::currVertexIndex;
    ImGui::Text("Singularity: ");
    ImGui::SameLine(); // put label on the left
    ImGui::InputDouble("", idx == -1 ? &INDEX : &RHO[idx], -std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity(), "%0.2f");

    if (ImGui::Button("Compute <phi>")) {
        Vector<double> phi = connections.computeConnections(RHO);
        if (phi.norm() != 0) {
            updateDirectionField(phi, true);
            redraw();
        }
    }
    if (ImGui::Button("Reset")) {
        polyscope::state::subset.vertices.clear();
        polyscope::state::currVertexIndex = -1;
        RHO = Vector<double>::Zero(mesh->nVertices());
        directionField->setEnabled(false);
        redraw();
    }
}


int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("15-458 HW6");
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

    // Load mesh
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();

    EULER_CHARACTERISTIC = geometry->eulerCharacteristic();

    if (mesh->nBoundaryLoops() > 0) {
        std::cerr << "Please load a mesh without boundary." << std::endl;
        return EXIT_FAILURE;
    }

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = functionCallback;

    // Add mesh to GUI
    psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(filepath), geometry->inputVertexPositions,
                                            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

    // Initialize operators.
    flipZ();

    // Add visualization options.
    psMesh->setSmoothShade(true);
    psMesh->setSurfaceColor({1.0, 0.45, 0.0}); // orange
    lengthScale = geometry->meanEdgeLength();
    vertexRadius = lengthScale * 0.2;

    std::srand(std::time(0));
    connections = TrivialConnections(mesh, geometry);
    RHO = Vector<double>::Zero(mesh->nVertices());
    // initialize to something
    currVert =
        psMesh->addSurfaceGraphQuantity("current vertex", std::vector<Vector3>(), std::vector<std::array<size_t, 2>>());
    updateDirectionField(Vector<double>::Zero(mesh->nEdges()), false);
    // Get indices for element picking
    polyscope::state::facePickIndStart = mesh->nVertices();
    polyscope::state::edgePickIndStart = polyscope::state::facePickIndStart + mesh->nFaces();
    polyscope::state::halfedgePickIndStart = polyscope::state::edgePickIndStart + mesh->nEdges();


    // Give control to the polyscope gui
    polyscope::show();

    delete mesh;
    delete geometry;

    return EXIT_SUCCESS;
}
