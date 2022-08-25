// DISCRETE CURVATURES AND NORMALS

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "colormap.h"
#include <map>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

// Some global variables
double TOTAL_ANGLE_DEFECT;
size_t EULER_CHARACTERISTIC;
// normal vectors
enum normalType { nOFF, nEQUAL, nTIP, nSPHERE, nAREA, nMEAN, nGAUSS };
std::vector<std::string> nNames = {"Off",
                                   "Equally weighted",
                                   "Tip angle weighted",
                                   "Sphere inscribed",
                                   "Area weighted (AN)",
                                   "Mean curvature (HN)",
                                   "Gauss curvature (KN)"};
std::map<int, std::vector<Vector3>> nVectors;
polyscope::SurfaceVectorQuantity* normalVectors;
std::vector<Vector3> EW_vectors;
std::vector<Vector3> TW_vectors;
std::vector<Vector3> SI_vectors;
std::vector<Vector3> AN_vectors;
std::vector<Vector3> HN_vectors;
std::vector<Vector3> KN_vectors;
// vertex colors
enum shadeType { sSHADED, sAREA, sMEAN, sGAUSS, sKMAX, sKMIN };
std::vector<std::string> sNames = {"Off", "A", "H", "K", "k max", "k min"};
std::map<int, std::vector<std::array<double, 3>>> sColors;
polyscope::SurfaceVertexColorQuantity* vertexColors;
std::vector<std::array<double, 3>> shaded_colors;
std::vector<std::array<double, 3>> A_colors;
std::vector<std::array<double, 3>> H_colors;
std::vector<std::array<double, 3>> K_colors;
std::vector<std::array<double, 3>> kmax_colors;
std::vector<std::array<double, 3>> kmin_colors;


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

/* Called once at start of program. */
void computeNormals() {

    for (Vertex v : mesh->vertices()) {
        EW_vectors.push_back(geometry->vertexNormalEquallyWeighted(v));
        TW_vectors.push_back(geometry->vertexNormalAngleWeighted(v));
        SI_vectors.push_back(geometry->vertexNormalSphereInscribed(v));
        AN_vectors.push_back(geometry->vertexNormalAreaWeighted(v));
        HN_vectors.push_back(geometry->vertexNormalMeanCurvature(v));
        KN_vectors.push_back(geometry->vertexNormalGaussianCurvature(v));
    }
    nVectors[nEQUAL] = EW_vectors;
    nVectors[nTIP] = TW_vectors;
    nVectors[nSPHERE] = SI_vectors;
    nVectors[nAREA] = AN_vectors;
    nVectors[nMEAN] = HN_vectors;
    nVectors[nGAUSS] = KN_vectors;
}

/* Computed once at start of program. */
void computeShaded() {

    // Get values
    Vector<double> A(mesh->nVertices());
    double A_max = 0;
    Vector<double> H(mesh->nVertices());
    double H_max = 0;
    Vector<double> K(mesh->nVertices());
    double K_max = 0;
    Vector<double> kmin(mesh->nVertices());
    double kmin_max = 0;
    Vector<double> kmax(mesh->nVertices());
    double kmax_max = 0;
    for (Vertex v : mesh->vertices()) {
        size_t i = v.getIndex();
        A[i] = geometry->circumcentricDualArea(v);
        H[i] = geometry->scalarMeanCurvature(v);
        K[i] = geometry->angleDefect(v);
        double area = geometry->barycentricDualArea(v);
        std::pair<double, double> pc = geometry->principalCurvatures(v);
        double k1 = std::min(pc.first, pc.second);
        double k2 = std::max(pc.first, pc.second);
        kmin[i] = k1 * area;
        kmax[i] = k2 * area;
        // determine color map range
        A_max = std::max(A[i], A_max);
        H_max = std::max(abs(H[i]), H_max);
        K_max = std::max(abs(K[i]), K_max);
        kmin_max = std::max(abs(kmin[i]), kmin_max);
        kmax_max = std::max(abs(kmax[i]), kmax_max);
    }
    // Before implementing function, have mesh just be a neutral color (white)
    if (A_max == 0) {
        A_max = 1e-3;
    }
    if (H_max == 0) {
        H_max = 1e-3;
    }
    if (K_max == 0) {
        K_max = 1e-3;
    }
    if (kmin_max == 0) {
        kmin_max = 1e-3;
    }
    if (kmax_max == 0) {
        kmax_max = 1e-3;
    }
    // Determine colors
    for (size_t i = 0; i < mesh->nVertices(); i++) {
        shaded_colors.push_back({1.0, 0.45, 0.0}); // orange
        A_colors.push_back(mapToColor(A[i], -A_max, A_max, "seismic"));
        H_colors.push_back(mapToColor(H[i], -H_max, H_max, "seismic"));
        K_colors.push_back(mapToColor(K[i], -K_max, K_max, "seismic"));
        kmin_colors.push_back(mapToColor(kmin[i], -kmin_max, kmin_max, "seismic"));
        kmax_colors.push_back(mapToColor(kmax[i], -kmax_max, kmax_max, "seismic"));
    }
    sColors[sSHADED] = shaded_colors;
    sColors[sAREA] = A_colors;
    sColors[sMEAN] = H_colors;
    sColors[sGAUSS] = K_colors;
    sColors[sKMAX] = kmax_colors;
    sColors[sKMIN] = kmin_colors;
}

static const ImVec4 pressColor = ImColor::HSV(1. / 7.0f, 0.6f, 0.6f); // gold
static const ImVec4 releaseColor{0.35f, 0.61f, 0.49f, 0.62f};         // default green
static ImVec4 off_currColor = pressColor;
static ImVec4 EW_currColor = releaseColor;
static ImVec4 TA_currColor = releaseColor;
static ImVec4 SI_currColor = releaseColor;
static ImVec4 AN_currColor = releaseColor;
static ImVec4 HN_currColor = releaseColor;
static ImVec4 KN_currColor = releaseColor;
std::vector<ImVec4*> nState = {&off_currColor, &EW_currColor, &TA_currColor, &SI_currColor,
                               &AN_currColor,  &HN_currColor, &KN_currColor};

static ImVec4 shaded_currColor = pressColor;
static ImVec4 A_currColor = releaseColor;
static ImVec4 H_currColor = releaseColor;
static ImVec4 K_currColor = releaseColor;
static ImVec4 kmin_currColor = releaseColor;
static ImVec4 kmax_currColor = releaseColor;
std::vector<ImVec4*> sState = {&shaded_currColor, &A_currColor,    &H_currColor,
                               &K_currColor,      &kmin_currColor, &kmax_currColor};

void functionCallback() {

    ImGui::Text("Total angle defect: %0.1fpi", TOTAL_ANGLE_DEFECT / M_PI);
    ImGui::Text("Euler characteristic: %zu", EULER_CHARACTERISTIC);
    ImGui::Text("");

    ImGui::Text("Normals:");

    for (size_t i = 0; i < nNames.size(); i++) {
        ImGui::PushStyleColor(ImGuiCol_Button, *nState[i]);
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, *nState[i]);
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, *nState[i]);
        if (ImGui::Button(nNames[i].c_str())) {
            if (i != nOFF) {
                normalVectors = psMesh->addVertexVectorQuantity("Normals", nVectors[i]);
                normalVectors->setEnabled(true);
            } else {
                normalVectors->setEnabled(false);
            }
            for (size_t j = 0; j < nVectors.size(); j++) {
                *nState[j] = i == j ? pressColor : releaseColor;
            }
        }
        ImGui::PopStyleColor(3);
        if (i % 3 != 2 && i != nNames.size() - 1) {
            ImGui::SameLine();
        }
    }
    normalVectors->setVectorRadius(0.001);
    normalVectors->setVectorColor({0.0, 0.0, 1.0});

    ImGui::Text("Plot:");

    for (size_t i = 0; i < sColors.size(); i++) {
        ImGui::PushStyleColor(ImGuiCol_Button, *sState[i]);
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, *sState[i]);
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, *sState[i]);
        if (ImGui::Button(sNames[i].c_str())) {
            vertexColors = psMesh->addVertexColorQuantity("Plot", sColors[i]);
            for (size_t j = 0; j < sColors.size(); j++) {
                *sState[j] = i == j ? pressColor : releaseColor;
            }
        }
        ImGui::PopStyleColor(3);
        ImGui::SameLine();
    }
}


int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("15-458 HW2");
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
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(filepath);

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = functionCallback;

    // Add mesh to GUI
    psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(filepath), geometry->inputVertexPositions,
                                            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

    // Add visualization options.
    flipZ();
    psMesh->setSmoothShade(true);
    psMesh->setSurfaceColor({1.0, 0.45, 0.0});
    computeNormals();
    normalVectors = psMesh->addVertexVectorQuantity("Normals", EW_vectors);
    normalVectors->setEnabled(false);
    computeShaded();
    vertexColors = psMesh->addVertexColorQuantity("Plot", shaded_colors);
    vertexColors->setEnabled(true);

    // Initialize quantities.
    TOTAL_ANGLE_DEFECT = geometry->totalAngleDefect();
    EULER_CHARACTERISTIC = geometry->eulerCharacteristic();

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}