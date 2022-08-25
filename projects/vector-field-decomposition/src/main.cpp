// VECTOR FIELD DECOMPOSITION

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "harmonic-bases.h"
#include "hodge-decomposition.h"
#include "scalar-poisson-problem.h"
#include "tree-cotree.h"

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;
// so we can more easily pass these to different classes
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

// Some global variables
TreeCotree treeCotree;
HodgeDecomposition hodgeDecomposition;
HarmonicBases harmonicBases;
ScalarPoissonProblem SPP;
Vector<double> omega, dAlpha, deltaBeta, myGamma;
std::vector<Vector<double>> bases;

polyscope::SurfaceGraphQuantity* treeGraph;
polyscope::SurfaceGraphQuantity* cotreeGraph;
polyscope::SurfaceGraphQuantity* generatorsGraph;
polyscope::SurfaceVectorQuantity* omegaVis;
polyscope::SurfaceVectorQuantity* dAlphaVis;
polyscope::SurfaceVectorQuantity* deltaBetaVis;
polyscope::SurfaceVectorQuantity* gammaVis;
std::vector<polyscope::SurfaceVectorQuantity*> basesVis;

int currentBasis = 0;
size_t EULER_CHARACTERISTIC;
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

void showSelected() {
    // pass
}

size_t computeEulerCharacteristic() {
    return mesh->nVertices() - mesh->nEdges() + mesh->nFaces();
}

Vector3 centroid(Face f) {

    Halfedge he = f.halfedge();
    Vector3 a = geometry->inputVertexPositions[he.vertex()];
    Vector3 b = geometry->inputVertexPositions[he.next().vertex()];
    Vector3 c = geometry->inputVertexPositions[he.next().next().vertex()];
    return (a + b + c) / 3.0;
}


Vector<double> generateRandomOneForm() {

    // Generate random vector and scalar potentials.
    int V = mesh->nVertices();
    int n = std::max(2, (int)(V / 1000));
    Vector<double> rho1 = Vector<double>::Zero(V);
    Vector<double> rho2 = Vector<double>::Zero(V);
    for (int i = 0; i < n; i++) {
        rho1[rand() % V] = rand() % 1000 - 500;
        rho2[rand() % V] = rand() % 1000 - 500;
    }
    Vector<double> scalarPotential = SPP.solve(rho1);
    Vector<double> vectorPotential = SPP.solve(rho2);

    // Compute per-face field
    std::map<Face, Vector3> field;
    for (Face f : mesh->faces()) {
        double A = geometry->faceArea(f);
        Vector3 N = geometry->faceNormal(f);
        Vector3 C = centroid(f);
        field[f] = {0, 0, 0};

        // Add exact and coexact components
        for (Halfedge he : f.adjacentHalfedges()) {
            size_t i = he.next().tipVertex().getIndex();
            Vector3 e = geometry->halfedgeVector(he);
            Vector3 eT = cross(N, e);
            field[f] += eT * scalarPotential[i] / (2 * A);
            field[f] += e * vectorPotential[i] / (2 * A);
        }

        // Add harmonic component.
        Vector3 u = {-C[2], 0.0, C[0]};
        u -= N * dot(u, N);
        u = u.normalize();
        field[f] += u * 0.3;
    }

    // Compute 1-form
    Vector<double> form(mesh->nEdges());
    for (Edge e : mesh->edges()) {
        Halfedge he = e.halfedge(); // should always be an interior halfedge
        Vector3 f1 = field[he.face()];
        Vector3 f2 = he.twin().isInterior() ? field[he.twin().face()] : Vector3{0, 0, 0};
        form[e.getIndex()] = 0.5 * dot(f1 + f2, geometry->halfedgeVector(he));
    }
    return form;
}


/*
 * Create surface graph quantity representing the tree-cotree.
 */
void buildTreeCotreeGraph() {

    // Build primal spanning tree.
    std::vector<Vector3> vertPrimal;
    std::vector<std::array<size_t, 2>> indPrimal;
    if (treeCotree.vertexParent.size() > 0) {
        for (Vertex v : mesh->vertices()) {
            vertPrimal.push_back(geometry->inputVertexPositions[v]);
            indPrimal.push_back({v.getIndex(), treeCotree.vertexParent[v].getIndex()});
        }
    }

    // Build dual spanning tree.
    std::vector<Vector3> vertDual;
    std::vector<std::array<size_t, 2>> indDual;
    if (treeCotree.faceParent.size() > 0) {
        for (Face f : mesh->faces()) {
            vertDual.push_back(centroid(f));
            indDual.push_back({f.getIndex(), treeCotree.faceParent[f].getIndex()});
        }
    }

    treeGraph = psMesh->addSurfaceGraphQuantity("Tree", vertPrimal, indPrimal);
    cotreeGraph = psMesh->addSurfaceGraphQuantity("Cotree", vertDual, indDual);
    treeGraph->setColor({0, 0, 1});
    cotreeGraph->setColor({0, 1, 0});
    treeGraph->setRadius(lengthScale * 0.04);
    cotreeGraph->setRadius(lengthScale * 0.04);
    treeGraph->setEnabled(false);
    cotreeGraph->setEnabled(false);
}

/*
 * Create surface graph quantity representing the generators.
 */
void buildGeneratorsGraph() {

    std::vector<Vector3> verts;
    std::vector<std::array<size_t, 2>> inds;
    size_t offset = 0;
    for (size_t i = 0; i < treeCotree.generators.size(); i++) {
        size_t M = treeCotree.generators[i].size();
        // Add interior vertices and edges
        for (size_t j = 0; j < M; j++) {
            Halfedge he = treeCotree.generators[i][j];
            Vector3 a = geometry->inputVertexPositions[he.tailVertex()];
            Vector3 b = geometry->inputVertexPositions[he.tipVertex()];
            verts.push_back(centroid(he.face()));
            verts.push_back((a + b) * 0.5);
            verts.push_back(centroid(he.twin().face()));
            inds.push_back({offset + 3 * j, offset + 3 * j + 1});
            inds.push_back({offset + 3 * j + 1, offset + 3 * j + 2});
        }
        offset += 3 * M;
    }
    generatorsGraph = psMesh->addSurfaceGraphQuantity("Generators", verts, inds);
    generatorsGraph->setColor({1, 0, 0});
    generatorsGraph->setEnabled(false);
}

/*
 * Add visualizations for harmonic bases.
 */
void buildBasesVis() {

    if (bases.size() == 0) {
        basesVis.push_back(psMesh->addOneFormIntrinsicVectorQuantity("Empty", Vector<double>::Zero(mesh->nEdges()),
                                                                     polyscopeEdgeOrientations(*mesh)));
        return;
    }

    for (size_t i = 0; i < bases.size(); i++) {
        std::ostringstream stringStream;
        stringStream << "gamma_" << i;
        std::string basisName = stringStream.str();
        polyscope::SurfaceVectorQuantity* gamma_i =
            psMesh->addOneFormIntrinsicVectorQuantity(basisName, bases[i], polyscopeEdgeOrientations(*mesh));
        gamma_i->setEnabled(false);
        gamma_i->setVectorColor({0, 0, 0});
        basesVis.push_back(gamma_i);
    }
}

/*
 * Update quantities when user generates new <omega>.
 */
void updateQuantities() {

    omega = generateRandomOneForm();
    dAlpha = hodgeDecomposition.computeExactComponent(omega);
    deltaBeta = hodgeDecomposition.computeCoExactComponent(omega);
    myGamma = hodgeDecomposition.computeHarmonicComponent(omega, dAlpha, deltaBeta);
    if (myGamma.norm() / omega.norm() < 1e-6) {
        myGamma = Vector<double>::Zero(mesh->nEdges());
    }

    // So program doesn't give error if tasks haven't been implemented yet
    if (dAlpha.size() != (int)mesh->nEdges()) {
        dAlpha = Vector<double>::Zero(mesh->nEdges());
    }
    if (deltaBeta.size() != (int)mesh->nEdges()) {
        deltaBeta = Vector<double>::Zero(mesh->nEdges());
    }
    if (myGamma.size() != (int)mesh->nEdges()) {
        myGamma = Vector<double>::Zero(mesh->nEdges());
    }

    // Update visualizations.
    omegaVis = psMesh->addOneFormIntrinsicVectorQuantity("omega", omega, polyscopeEdgeOrientations(*mesh));
    dAlphaVis = psMesh->addOneFormIntrinsicVectorQuantity("dAlpha", dAlpha, polyscopeEdgeOrientations(*mesh));
    deltaBetaVis = psMesh->addOneFormIntrinsicVectorQuantity("deltaBeta", deltaBeta, polyscopeEdgeOrientations(*mesh));
    gammaVis = psMesh->addOneFormIntrinsicVectorQuantity("gamma", myGamma, polyscopeEdgeOrientations(*mesh));

    omegaVis->setVectorColor({0, 0, 0});
    dAlphaVis->setVectorColor({0, 0, 0});
    deltaBetaVis->setVectorColor({0, 0, 0});
    gammaVis->setVectorColor({0, 0, 0});
}


// have not yet figured out how to get characters to render properly in ImGui
enum quantityName { OFF, OMEGA, D_ALPHA, DELTA_BETA, GAMMA, TREE, GENERATORS, BASES };
const char* items[] = {
    "Off",         "1-form <omega>", "Exact part d<alpha>",     "Coexact part delta<beta>", "Harmonic part <gamma>",
    "Tree cotree", "Generators",     "Harmonic bases <gamma>_i"};
static const char* current_item = "1-form <omega>";

/*
 * Enable/disable visualizations depending on which one is selected.
 * This would be a lot more efficient if all the vis structures were in an array
 */
void updateVis(int activate = 0) {

    switch (activate) {
    case OFF:
        treeGraph->setEnabled(false);
        cotreeGraph->setEnabled(false);
        generatorsGraph->setEnabled(false);
        omegaVis->setEnabled(false);
        dAlphaVis->setEnabled(false);
        deltaBetaVis->setEnabled(false);
        gammaVis->setEnabled(false);
        basesVis[currentBasis]->setEnabled(false);
        break;
    case OMEGA:
        treeGraph->setEnabled(false);
        cotreeGraph->setEnabled(false);
        generatorsGraph->setEnabled(false);
        omegaVis->setEnabled(true);
        dAlphaVis->setEnabled(false);
        deltaBetaVis->setEnabled(false);
        gammaVis->setEnabled(false);
        basesVis[currentBasis]->setEnabled(false);
        break;
    case D_ALPHA:
        treeGraph->setEnabled(false);
        cotreeGraph->setEnabled(false);
        generatorsGraph->setEnabled(false);
        omegaVis->setEnabled(false);
        dAlphaVis->setEnabled(true);
        deltaBetaVis->setEnabled(false);
        gammaVis->setEnabled(false);
        basesVis[currentBasis]->setEnabled(false);
        break;
    case DELTA_BETA:
        treeGraph->setEnabled(false);
        cotreeGraph->setEnabled(false);
        generatorsGraph->setEnabled(false);
        omegaVis->setEnabled(false);
        dAlphaVis->setEnabled(false);
        deltaBetaVis->setEnabled(true);
        gammaVis->setEnabled(false);
        basesVis[currentBasis]->setEnabled(false);
        break;
    case GAMMA:
        treeGraph->setEnabled(false);
        cotreeGraph->setEnabled(false);
        generatorsGraph->setEnabled(false);
        omegaVis->setEnabled(false);
        dAlphaVis->setEnabled(false);
        deltaBetaVis->setEnabled(false);
        gammaVis->setEnabled(true);
        basesVis[currentBasis]->setEnabled(false);
        break;
    case TREE:
        treeGraph->setEnabled(true);
        cotreeGraph->setEnabled(true);
        generatorsGraph->setEnabled(false);
        omegaVis->setEnabled(false);
        dAlphaVis->setEnabled(false);
        deltaBetaVis->setEnabled(false);
        gammaVis->setEnabled(false);
        basesVis[currentBasis]->setEnabled(false);
        break;
    case GENERATORS:
        treeGraph->setEnabled(false);
        cotreeGraph->setEnabled(false);
        generatorsGraph->setEnabled(true);
        omegaVis->setEnabled(false);
        dAlphaVis->setEnabled(false);
        deltaBetaVis->setEnabled(false);
        gammaVis->setEnabled(false);
        basesVis[currentBasis]->setEnabled(false);
        break;
    case BASES:
        treeGraph->setEnabled(false);
        cotreeGraph->setEnabled(false);
        generatorsGraph->setEnabled(false);
        omegaVis->setEnabled(false);
        dAlphaVis->setEnabled(false);
        deltaBetaVis->setEnabled(false);
        gammaVis->setEnabled(false);
        basesVis[currentBasis]->setEnabled(true);
        break;
    }
}

/*
 * User-defined buttons
 */
void functionCallback() {

    ImGui::Text("Euler characteristic: %zu", EULER_CHARACTERISTIC);
    ImGui::Text("");

    if (ImGui::Button("Generate new <omega>")) {
        updateQuantities();
        current_item = "1-form <omega>";
        updateVis(OMEGA);
        polyscope::requestRedraw();
    }

    if (ImGui::Button("Cycle to next <gamma>_i")) {
        if (bases.size() > 0) {
            basesVis[currentBasis]->setEnabled(false);
            currentBasis = (currentBasis + 1) % bases.size();
            current_item = "Harmonic bases <gamma>_i";
            updateVis(BASES);
            polyscope::requestRedraw();
        }
    }
    ImGui::Text("");

    ImGui::Text("Plot:");
    if (ImGui::BeginCombo("##plot", current_item)) {
        ImGui::SetItemDefaultFocus();
        for (int n = 0; n < IM_ARRAYSIZE(items); n++) {
            bool is_selected = (current_item == items[n]);
            if (ImGui::Selectable(items[n], is_selected)) {
                updateVis(n);
                current_item = items[n];
            }
            if (is_selected) {
                ImGui::SetItemDefaultFocus();
            }
        }
        polyscope::requestRedraw();
    }
    ImGui::EndCombo(); // dropdown menu only works if EndCombo() is outside if; I'm not sure why.
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

    EULER_CHARACTERISTIC = computeEulerCharacteristic();

    std::srand(std::time(0));
    treeCotree = TreeCotree(mesh, geometry);
    treeCotree.buildGenerators();
    hodgeDecomposition = HodgeDecomposition(mesh, geometry);
    harmonicBases = HarmonicBases(mesh, geometry);
    bases = harmonicBases.compute(treeCotree.generators, hodgeDecomposition);
    buildTreeCotreeGraph(); // visualize trees
    buildGeneratorsGraph(); // visualize generators
    buildBasesVis();        // visualize bases
    SPP = ScalarPoissonProblem(mesh, geometry);
    updateQuantities(); // generate random 1-form
    omegaVis->setEnabled(true);

    // Give control to the polyscope gui
    polyscope::show();

    delete mesh;
    delete geometry;

    return EXIT_SUCCESS;
}
