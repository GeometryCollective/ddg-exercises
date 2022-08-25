// GEOMETRIC FLOW

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "mean-curvature-flow.h"
#include "modified-mean-curvature-flow.h"

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

// Some global variables
float TIMESTEP = 0.001;
VertexData<Vector3> ORIG_VPOS; // original vertex positions
Vector3 CoM;                   // original center of mass
MeanCurvatureFlow MCF;
ModifiedMeanCurvatureFlow ModMCF;


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

void showSelected() {
    // pass
}

void redraw() {
    psMesh->updateVertexPositions(geometry->inputVertexPositions);
    polyscope::requestRedraw();
}

/*
 * User-defined buttons
 */
void functionCallback() {

    if (ImGui::Button("Mean curvature flow")) {
        MCF.integrate(TIMESTEP);
        geometry->normalize(CoM, false);
        redraw();
    }
    if (ImGui::Button("Modified mean curvature flow")) {
        ModMCF.integrate(TIMESTEP);
        geometry->normalize(CoM, false);
        redraw();
    }
    if (ImGui::Button("Reset")) {
        geometry->inputVertexPositions = ORIG_VPOS;
        psMesh->updateVertexPositions(ORIG_VPOS);
        polyscope::requestRedraw();
    }
    ImGui::SliderFloat("Timestep", &TIMESTEP, 0.001, 0.1);
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

    // Load mesh
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = functionCallback;

    // Add mesh to GUI
    psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(filepath), geometry->inputVertexPositions,
                                            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

    // Initialize operators.
    flipZ();
    ORIG_VPOS = geometry->inputVertexPositions;
    CoM = geometry->centerOfMass();
    MCF = MeanCurvatureFlow(mesh, geometry);
    ModMCF = ModifiedMeanCurvatureFlow(mesh, geometry);

    // Add visualization options.
    psMesh->setSmoothShade(true);
    psMesh->setSurfaceColor({1.0, 0.45, 0.0}); // orange

    // Give control to the polyscope gui
    polyscope::show();

    delete mesh;
    delete geometry;

    return EXIT_SUCCESS;
}