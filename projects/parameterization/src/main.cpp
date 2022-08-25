// PARAMETERIZATION

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "distortion.h"
#include "spectral-conformal-parameterization.h"

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
polyscope::SurfaceParameterizationQuantity* checkerboard;
Vector3 CoM;                        // original center of mass, for re-centering purposes
VertexData<Vector3> ORIGINAL;       // original vertex positions
VertexData<Vector3> SCP_MESH;       // vertex positions of SCP
VertexData<Vector2> SCP_FLATTENING; // SCP solution
// colors-per-vertex to show QC error of SCP
std::vector<std::array<double, 3>> SCP_QC_COLORS;
// colors-per-vertex to show area scaling of SCP
std::vector<std::array<double, 3>> SCP_AREA_COLORS;
bool DISPLAY_FLAT = false;


void flipZ() {
    // Don't know how best way to change <polyscope::view::viewMat>/homeview
    // externally Rotate mesh 180 deg about up-axis on startup
    glm::mat4x4 rot = glm::rotate(glm::mat4x4(1.0f), static_cast<float>(PI), glm::vec3(0, 1, 0));
    for (Vertex v : mesh->vertices()) {
        Vector3 vec = geometry->inputVertexPositions[v];
        glm::vec4 rvec = {vec[0], vec[1], vec[2], 1.0};
        rvec = rot * rvec;
        geometry->inputVertexPositions[v] = {rvec[0], rvec[1], rvec[2]};
    }
    psMesh->updateVertexPositions(geometry->inputVertexPositions);
}

void showSelected() {
    // pass
}

void redraw() {
    psMesh->updateVertexPositions(geometry->inputVertexPositions);
    polyscope::requestRedraw();
}

/*
 * Convert map to mesh positions.
 */
VertexData<Vector3> mapToMeshData(VertexData<Vector2>& flattening) {

    VertexData<Vector3> flat = geometry->inputVertexPositions;
    for (Vertex v : mesh->vertices()) {
        flat[v] = Vector3{flattening[v][0], flattening[v][1], 0.0};
    }
    return flat;
}

static const ImVec4 pressColor = ImColor::HSV(1. / 7.0f, 0.6f, 0.6f); // gold
static const ImVec4 releaseColor{0.35f, 0.61f, 0.49f, 0.62f};         // default green
static ImVec4 shaded_currColor = pressColor;
static ImVec4 qc_currColor = releaseColor;
static ImVec4 area_currColor = releaseColor;

void functionCallback() {

    ImGui::Text("Plot: ");
    ImGui::SameLine();

    ImGui::PushStyleColor(ImGuiCol_Button, shaded_currColor);
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, shaded_currColor);
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, shaded_currColor);
    if (ImGui::Button("Shaded")) {
        psMesh->setSurfaceColor({1.0, 0.45, 0.0}); // orange
        checkerboard->setEnabled(true);

        shaded_currColor = pressColor;
        qc_currColor = releaseColor;
        area_currColor = releaseColor;
    }
    ImGui::PopStyleColor(3);
    ImGui::SameLine();

    ImGui::PushStyleColor(ImGuiCol_Button, qc_currColor);
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, qc_currColor);
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, qc_currColor);
    if (ImGui::Button("Conformal error")) {
        polyscope::SurfaceColorQuantity* CE = psMesh->addVertexColorQuantity("Conformal error", SCP_QC_COLORS);
        CE->setEnabled(true);
        checkerboard->setEnabled(false);

        shaded_currColor = releaseColor;
        qc_currColor = pressColor;
        area_currColor = releaseColor;
    }
    ImGui::PopStyleColor(3);
    ImGui::SameLine();

    ImGui::PushStyleColor(ImGuiCol_Button, area_currColor);
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, area_currColor);
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, area_currColor);
    if (ImGui::Button("Area scaling")) {
        polyscope::SurfaceColorQuantity* AS = psMesh->addVertexColorQuantity("Area scaling", SCP_AREA_COLORS);
        AS->setEnabled(true);
        checkerboard->setEnabled(false);

        shaded_currColor = releaseColor;
        qc_currColor = releaseColor;
        area_currColor = pressColor;
    }
    ImGui::PopStyleColor(3);

    if (ImGui::Checkbox("Parameterization", &DISPLAY_FLAT)) {
        if (DISPLAY_FLAT) {
            // Display SCP parameterization
            geometry->inputVertexPositions = SCP_MESH;
            geometry->normalize(CoM, true);
            polyscope::view::flyToHomeView();
            polyscope::view::style = polyscope::view::NavigateStyle::Planar;
        } else {
            // Display original
            geometry->inputVertexPositions = ORIGINAL;
            polyscope::view::style = polyscope::view::NavigateStyle::Turntable;
        }
        redraw();
    }
}

/*
 * Add checkerboard pattern to mesh
 */
void addCheckerboard(VertexData<Vector2>& flattening) {

    std::vector<std::array<double, 2>> V_uv(mesh->nVertices());
    for (Vertex v : mesh->vertices()) {
        V_uv[v.getIndex()] = {flattening[v][0], flattening[v][1]};
    }
    checkerboard = psMesh->addVertexParameterizationQuantity("checkerboard", V_uv);
    checkerboard->setEnabled(true);
    checkerboard->setStyle(polyscope::ParamVizStyle::CHECKER);
    checkerboard->setCheckerColors(std::make_pair(glm::vec3{1.0, 0.45, 0.0}, glm::vec3{0.55, 0.27, 0.07}));
    checkerboard->setCheckerSize(0.002);
}

// Initalize objects for parameterization -- done once on program startup
void initializeParameterizations() {

    // Store initial mesh positions.
    ORIGINAL = geometry->inputVertexPositions;

    // Solve for parameterization and store mesh positions.
    SpectralConformalParameterization SCP = SpectralConformalParameterization(mesh, geometry);
    SCP_FLATTENING = SCP.flatten();
    SCP_MESH = mapToMeshData(SCP_FLATTENING);

    // Compute and store mesh colors for visualization.
    Distortion DIS = Distortion(mesh, geometry);
    DIS.computeQuasiConformalError(SCP_QC_COLORS, SCP_FLATTENING);
    DIS.computeAreaScaling(SCP_AREA_COLORS, SCP_FLATTENING);
}


int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("15-458 HW4");
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
    std::string filepath = "../../../input/face.obj";
    if (inputFilename) {
        filepath = args::get(inputFilename);
    }

    MESHNAME = polyscope::guessNiceNameFromPath(filepath);

    // Load mesh
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();

    if (mesh->nBoundaryLoops() == 0) {
        std::cerr << "Please load a mesh with at least 1 boundary loop" << std::endl;
        return EXIT_SUCCESS;
    }

    // Visualization stuff
    CoM = geometry->centerOfMass();
    geometry->normalize(CoM, true); // rescale to unit radius

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = functionCallback;

    // Add mesh to GUI
    psMesh = polyscope::registerSurfaceMesh(MESHNAME, geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));

    // Set up mesh for viewing
    flipZ();
    psMesh->setSmoothShade(true);
    psMesh->setSurfaceColor({1.0, 0.45, 0.0}); // orange

    initializeParameterizations();
    addCheckerboard(SCP_FLATTENING);

    // Give control to the polyscope gui
    polyscope::show();

    delete mesh;
    delete geometry;

    return EXIT_SUCCESS;
}