// DISCRETE EXTERIOR CALCULUS

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "colormap.h"
#include "solvers.h"
#include <algorithm>
#include <map>

// ** taken from Polyscope gl_engine.h -- not sure how to avoid using OpenGL for texture loading
#include "stb_image.h"
#ifdef __APPLE__
#define GLFW_INCLUDE_GLCOREARB
#include "GLFW/glfw3.h"
#else
#include "glad/glad.h"
// glad must come first
#include "GLFW/glfw3.h"
#endif

#ifdef _WIN32
#undef APIENTRY
#define GLFW_EXPOSE_NATIVE_WIN32
#define GLFW_EXPOSE_NATIVE_WGL
#include <GLFW/glfw3native.h>
#endif
// **

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

std::unique_ptr<ManifoldSurfaceMesh> mesh_dual;
std::unique_ptr<VertexPositionGeometry> geometry_dual;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;
polyscope::SurfaceMesh* psMesh_dual;

// Some global variables
enum formName { PRIMAL_0_FORM, PRIMAL_1_FORM, PRIMAL_2_FORM, DUAL_0_FORM, DUAL_1_FORM, DUAL_2_FORM };

std::map<int, std::string> formNameMap = {{PRIMAL_0_FORM, "Primal 0-form"}, {PRIMAL_1_FORM, "Primal 1-form"},
                                          {PRIMAL_2_FORM, "Primal 2-form"}, {DUAL_0_FORM, "Dual 0-form"},
                                          {DUAL_1_FORM, "Dual 1-form"},     {DUAL_2_FORM, "Dual 2-form"}};

struct DiagramImage {
    int width = 0, height = 0;
    GLuint texture = 0;
};

DiagramImage diagramImages[6] = {};

int currentFormName = PRIMAL_0_FORM;
std::map<Halfedge, size_t> heMap; // maps primal exterior halfedges to dual vertex indices
std::map<Vertex, size_t> vMap;    // maps primal boundary vertices to dual vertex indices
std::map<Edge, Edge> eMap;        // maps primal edges to dual edges
Vector<double> currentForm;
std::vector<std::array<double, 3>> formColor;
std::vector<char> EDGE_ORIENTATIONS;
std::vector<char> DUAL_EDGE_ORIENTATIONS;
polyscope::SurfaceVectorQuantity* formColorEdge;
polyscope::SurfaceVectorQuantity* formColorDualEdge;
double lengthScalePrimal, lengthScaleDual;

SparseMatrix<double> D0;
SparseMatrix<double> D1;
SparseMatrix<double> D0T;
SparseMatrix<double> D1T;
SparseMatrix<double> H0;
SparseMatrix<double> H1;
SparseMatrix<double> H2;
SparseMatrix<double> H0inv;
SparseMatrix<double> H1inv;
SparseMatrix<double> H2inv;


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

/*
 * Compute the midpoint of a halfedge.
 */
Vector3 midpoint(Halfedge he) {

    Vector3 p1 = geometry->inputVertexPositions[he.tailVertex()];
    Vector3 p2 = geometry->inputVertexPositions[he.tipVertex()];
    return (p1 + p2) / 2;
}

/*
 * Computes the circumcenter of a face.
 */
Vector3 circumcenter(Face f) {

    assert(!f.isBoundaryLoop());
    // Get the 3 vertices of the face
    Halfedge he = f.halfedge();
    Vector3 a = geometry->inputVertexPositions[he.vertex()];
    Vector3 b = geometry->inputVertexPositions[he.next().vertex()];
    Vector3 c = geometry->inputVertexPositions[he.next().next().vertex()];

    Vector3 ac = c - a;
    Vector3 ab = b - a;
    Vector3 w = cross(ab, ac);

    Vector3 u = cross(w, ab) * ac.norm() * ac.norm();
    Vector3 v = cross(ac, w) * ab.norm() * ab.norm();
    Vector3 x = (u + v) * (1.0 / (2.0 * w.norm() * w.norm()));
    return x + a;
}

/*
 * Create dual mesh.
 */
void createDualMesh() {

    // map halfedges in original mesh to vertex indices in dual mesh
    size_t addl = 1;
    for (Halfedge he : mesh->exteriorHalfedges()) {
        heMap[he] = mesh->nFaces() + addl - 1;
        addl += 1;
    }
    // map primal boundary vertices to dual vertices
    addl = mesh->nFaces() + heMap.size();
    for (Halfedge he : mesh->exteriorHalfedges()) {
        vMap[he.vertex()] = addl;
        addl += 1;
    }

    // a list of 0-indexed vertices incident on each face in CCW order
    std::vector<std::vector<size_t>> polygons;
    for (Vertex v : mesh->vertices()) {
        std::vector<size_t> poly;
        if (v.isBoundary()) {
            for (Halfedge he : v.outgoingHalfedges()) {
                // get midpoint of one boundary halfedge
                if (!he.twin().isInterior()) {
                    // I thought these should be switched in order... but this seems to work
                    poly.push_back(he.face().getIndex());
                    poly.push_back(heMap[he.twin()]);
                    poly.push_back(vMap[v]);
                }
                // get midpoint of other boundary halfedge
                else if (!he.isInterior()) {
                    poly.push_back(heMap[he]);
                } else {
                    poly.push_back(he.face().getIndex());
                }
            }
        } else {
            for (Face f : v.adjacentFaces()) {
                poly.push_back(f.getIndex());
            }
        }
        polygons.push_back(poly);
    }
    mesh_dual.reset(new ManifoldSurfaceMesh(polygons));

    // Get vertex positions
    std::unique_ptr<VertexPositionGeometry> dual_geo(new VertexPositionGeometry(*mesh_dual));
    for (Face f : mesh->faces()) {
        dual_geo->inputVertexPositions[f.getIndex()] = circumcenter(f);
    }
    for (Halfedge he : mesh->exteriorHalfedges()) {
        // An extra vertex in the middle of each exterior halfedge
        dual_geo->inputVertexPositions[heMap[he]] = midpoint(he);
        // and an extra vertex for every boundary vertex
        dual_geo->inputVertexPositions[vMap[he.vertex()]] = geometry->inputVertexPositions[he.vertex()];
    }
    geometry_dual = std::move(dual_geo);
}

// https://github.com/ocornut/imgui/wiki/Image-Loading-and-Displaying-Examples
bool LoadTextureFromFile(const char* filename, GLuint* out_texture, int* out_width, int* out_height) {
    // Load from file
    int image_width = 0;
    int image_height = 0;
    unsigned char* image_data = stbi_load(filename, &image_width, &image_height, NULL, 4);
    if (image_data == NULL) return false;

    // Create a OpenGL texture identifier
    GLuint image_texture;
    glGenTextures(1, &image_texture);
    glBindTexture(GL_TEXTURE_2D, image_texture);

    // Setup filtering parameters for display
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                    GL_CLAMP_TO_EDGE);                                   // This is required on WebGL for non
                                                                         // power-of-two textures
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE); // Same

// Upload pixels into texture
#if defined(GL_UNPACK_ROW_LENGTH) && !defined(__EMSCRIPTEN__)
    glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
#endif
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image_width, image_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image_data);
    stbi_image_free(image_data);

    *out_texture = image_texture;
    *out_width = image_width;
    *out_height = image_height;

    return true;
}

/*
 * Display the operator diagram.
 */
void showDiagram() {
    DiagramImage& image = diagramImages[currentFormName];

    if (image.texture == 0) {
        char filename[256];
        snprintf(filename, sizeof(filename), "../../../imgs/dec%d.png", currentFormName + 1);
        bool ret = LoadTextureFromFile(filename, &image.texture, &image.width, &image.height);
        IM_ASSERT(ret);
    }

    // inside function callback, ImGui::GetWindowSize() gets size of callback window
    float width = ImGui::GetWindowWidth();
    float imgwidth = width * 0.7;
    float imgheight = imgwidth * ((float)image.height / image.width);
    float windowHeight = imgheight + 40;
    float windowWidth = imgwidth + 40;
    static const ImVec4 transparent{1.0f, 1.0f, 1.0f, 0.0f};
    ImGui::PushStyleColor(ImGuiCol_WindowBg, transparent);
    ImGui::PushStyleColor(ImGuiCol_Border, transparent);
    ImGui::PushStyleColor(ImGuiCol_TitleBg, transparent);
    ImGui::PushStyleColor(ImGuiCol_TitleBgActive, transparent);
    ImGui::PushStyleColor(ImGuiCol_TitleBgCollapsed, transparent);
    ImGui::PushStyleColor(ImGuiCol_Text, transparent);
    ImGui::SetNextWindowPos(
        ImVec2(polyscope::view::windowWidth - windowWidth, polyscope::view::windowHeight - windowHeight));
    ImGui::SetNextWindowSize(ImVec2(windowWidth, windowHeight));
    if (ImGui::Begin("Diagram")) {
        ImGui::Image((void*)(intptr_t)image.texture, ImVec2(imgwidth, imgheight));
    }
    ImGui::End();
    ImGui::PopStyleColor(6);
}


double gaussian(double x, double a, double b) {
    return a * exp(-(x * x) / (b * b));
}

Vector<double> generateRandomFormOnVertices(bool scaleByArea = false) {

    // # of peaks is in [2, 10]
    int nPeaks = std::max(2, rand() % 10);
    // peaks are at random vertices
    std::vector<Vector3> c;
    for (int i = 0; i < nPeaks; i++) {
        c.push_back(geometry->inputVertexPositions[rand() % mesh->nVertices()]);
    }
    Vector<double> form = Vector<double>::Zero(mesh->nVertices());
    for (Vertex v : mesh->vertices()) {
        double sum = 0;
        for (int i = 0; i < nPeaks; i++) {
            sum += gaussian((c[i] - geometry->inputVertexPositions[v]).norm(), 1.0, 0.5);
        }
        if (scaleByArea) {
            sum *= geometry->barycentricDualArea(v);
        }
        form[v.getIndex()] = sum;
    }
    return form;
}

Vector<double> generateRandomFormOnFaces(bool scaleByArea = false) {

    int nPeaks = std::max(2, rand() % 10);
    std::vector<Vector3> c;
    for (int i = 0; i < nPeaks; i++) {
        c.push_back(circumcenter(mesh->face(rand() % mesh->nFaces())));
    }

    Vector<double> form = Vector<double>::Zero(mesh->nFaces());
    for (Face f : mesh->faces()) {
        double sum = 0;
        for (int i = 0; i < nPeaks; i++) {
            sum += gaussian((c[i] - circumcenter(f)).norm(), 1.0, 0.5);
        }
        if (scaleByArea) {
            sum *= geometry->faceArea(f);
        }
        form[f.getIndex()] = sum;
    }
    return form;
}

Vector<double> generateRandomFormOnEdges() {

    Vector<double> scalarPotential = generateRandomFormOnVertices();
    Vector<double> vectorPotential = generateRandomFormOnVertices();

    std::vector<Vector3> field;
    field.resize(mesh->nFaces());
    for (Face f : mesh->faces()) {
        field[f.getIndex()] = {0, 0, 0};
        double A = geometry->faceArea(f);
        Vector3 N = geometry->faceNormal(f);
        Vector3 C = circumcenter(f);

        for (Halfedge he : f.adjacentHalfedges()) {
            size_t i = he.next().tipVertex().getIndex();
            Vector3 e = geometry->halfedgeVector(he);
            Vector3 eT = cross(N, e);
            field[f.getIndex()] += eT * scalarPotential[i] / (2 * A);
            field[f.getIndex()] += e * vectorPotential[i] / (2 * A);
        }

        Vector3 u = {-C[1], 0.0, C[0]};
        u -= N * dot(u, N);
        field[f.getIndex()] += u;
    }

    Vector<double> form = Vector<double>::Zero(mesh->nEdges());
    for (Edge e : mesh->edges()) {
        Halfedge he = e.halfedge();
        Vector3 f1 = he.isInterior() ? field[he.face().getIndex()] : Vector3({0, 0, 0});
        Vector3 f2 = he.twin().isInterior() ? field[he.twin().face().getIndex()] : Vector3({0, 0, 0});
        Vector3 vec = geometry->halfedgeVector(he);
        form[e.getIndex()] = dot(f1 + f2, vec) * 0.5;
    }
    return form;
}

/*
 * Helper function for interpolateWachspressWhitney(). Computes the distance from a point <p> to the segment with
 * endpoint <a> and direction <u>.
 */
double pointToSegmentDistance(Vector3 p, Vector3 a, Vector3 u) {
    double t = std::max(0.0, std::min(1.0, dot(p - a, u) / u.norm2()));
    Vector3 v = a + t * u;
    return (p - v).norm();
}

/*
 * Instead of using Polyscope's addOneFormIntrinsicVectorQuantity(), for better control of the visualization.
 */
std::vector<Vector3> interpolateWhitney() {

    std::vector<Vector3> field(mesh->nFaces());
    for (Face f : mesh->faces()) {
        Halfedge he = f.halfedge();
        Vector3 pi = geometry->inputVertexPositions[he.vertex()];
        Vector3 pj = geometry->inputVertexPositions[he.next().vertex()];
        Vector3 pk = geometry->inputVertexPositions[he.next().next().vertex()];
        Vector3 eij = pj - pi;
        Vector3 ejk = pk - pj;
        Vector3 eki = pi - pk;
        double cij = currentForm[he.edge().getIndex()];
        double cjk = currentForm[he.next().edge().getIndex()];
        double cki = currentForm[he.next().next().edge().getIndex()];
        if (he.edge().halfedge() != he) {
            cij *= -1;
        }
        if (he.next().edge().halfedge() != he.next()) {
            cjk *= -1;
        }
        if (he.next().next().edge().halfedge() != he.next().next()) {
            cki *= -1;
        }
        double A = geometry->faceArea(f);
        Vector3 N = geometry->faceNormal(f);
        Vector3 a = (eki - ejk) * cij;
        Vector3 b = (eij - eki) * cjk;
        Vector3 c = (ejk - eij) * cki;
        Vector3 vec = cross(N, a + b + c) / (6 * A);
        double norm = (vec * lengthScalePrimal).norm();
        if (norm > lengthScalePrimal * 1.5) {
            vec *= lengthScalePrimal * 1.5 / norm;
        }
        field[f.getIndex()] = vec;
    }
    return field;
}
/*
 * Create dual 1-form visualization, since Polyscope's 1-form visualization function only takes triangles. Adapted from
 * JS version of this assignment.
 */
std::vector<Vector3> interpolateWachspressWhitney() {

    std::vector<Vector3> field(mesh->nVertices()); // == mesh_dual->nFaces()

    for (Vertex v : mesh->vertices()) {
        Vector3 p = geometry->inputVertexPositions[v];
        std::vector<Vector3> C;
        std::vector<double> w;

        if (v.isBoundary()) {
            continue;
        }
        for (Halfedge he : v.outgoingHalfedges()) {
            Vector3 f1 = circumcenter(he.twin().face());
            Vector3 f2 = circumcenter(he.face());
            Vector3 u = f2 - f1;
            double height = pointToSegmentDistance(p, f1, u);
            Vector3 normal = {u[1], -u[0], 0.0};
            C.push_back(normal / height);
            double wij = currentForm[he.edge().getIndex()];
            if (he.edge().halfedge() != he) {
                wij *= -1;
            }
            w.push_back(wij);
        }
        int n = C.size();
        for (int j = 0; j < n; j++) {
            int i = j == 0 ? n - 1 : j - 1; // handle mod of negative number
            int k = (j + 1) % n;
            field[v.getIndex()] += C[k] - w[j] * C[i];
        }
        Vector3 vec = field[v.getIndex()];
        double norm = (vec * lengthScaleDual).norm();
        if (norm > lengthScaleDual * 1.5) {
            vec *= lengthScaleDual * 1.5 / norm;
        }
        field[v.getIndex()] = vec;
    }
    return field;
}

/*
 * Determine the min and max values in a matrix.
 * Used for determining roughly what the range for the colormap should be.
 */
std::pair<double, double> minMaxFormValue(const Vector<double>& form) {

    double minVal = std::numeric_limits<double>::infinity();
    double maxVal = -std::numeric_limits<double>::infinity();

    if (currentFormName == PRIMAL_0_FORM || currentFormName == DUAL_2_FORM) {
        for (Vertex v : mesh->vertices()) {
            // don't bother with vertices on the boundary
            if (!v.isBoundary()) {
                double val = currentForm[v.getIndex()];
                if (currentFormName == DUAL_2_FORM) {
                    val /= geometry->barycentricDualArea(v);
                }
                minVal = std::min(val, minVal);
                maxVal = std::max(val, maxVal);
            }
        }
    } else if (currentFormName == PRIMAL_2_FORM || currentFormName == DUAL_0_FORM) {
        for (Face f : mesh->faces()) {
            bool hasBoundaryVertex = false;
            for (Vertex v : f.adjacentVertices()) {
                if (v.isBoundary()) {
                    hasBoundaryVertex = true;
                    break;
                }
            }
            // don't bother with faces on the boundary
            if (!hasBoundaryVertex) {
                double val = currentForm[f.getIndex()];
                if (currentFormName == PRIMAL_2_FORM) {
                    val /= geometry->faceArea(f);
                }
                minVal = std::min(val, minVal);
                maxVal = std::max(val, maxVal);
            }
        }
    } else {
        for (int i = 0; i < form.rows(); i++) {
            double val = currentForm[i];
            minVal = std::min(val, minVal);
            maxVal = std::max(val, maxVal);
        }
    }

    // so we don't get a non-range for the colormap
    double eps = 1e-5;
    if (abs(minVal) < eps && abs(maxVal) < eps) {
        for (int i = 0; i < currentForm.rows(); i++) {
            currentForm[i] = 0;
        }
    }
    if (abs(minVal) < eps) {
        minVal = 0;
    }
    if (abs(maxVal) < eps) {
        maxVal = eps;
    }

    return std::make_pair(minVal, maxVal);
}

void updateColorsPrimalForm() {

    // Determine the min and max values in the matrix.
    std::pair<double, double> minMaxForm = minMaxFormValue(currentForm);
    double minVal = minMaxForm.first;
    double maxVal = minMaxForm.second;
    formColor.clear();

    if (currentFormName == PRIMAL_0_FORM) {
        // Assign colors per vertex
        for (Vertex v : mesh->vertices()) {
            formColor.push_back(mapToColor(currentForm[v.getIndex()], minVal, maxVal, "coolwarm"));
        }
    } else if (currentFormName == PRIMAL_1_FORM) {
        // 1-forms visualized with vectors, so make background white
        for (Vertex v : mesh->vertices()) {
            formColor.push_back({1.0, 1.0, 1.0});
        }
        // Round very small values to 0 to avoid funky visualizations
        if (maxVal == 1e-5) {
            for (Edge e : mesh->edges()) {
                currentForm[e.getIndex()] = 0;
            }
        }
    } else if (currentFormName == PRIMAL_2_FORM) {
        // Assign colors per face
        for (Face f : mesh->faces()) {
            formColor.push_back(
                mapToColor(currentForm[f.getIndex()] / geometry->faceArea(f), minVal, maxVal, "coolwarm"));
        }
    }
}

void updateColorsDualForm() {

    // Determine the min and max values in the matrix.
    std::pair<double, double> minMaxForm = minMaxFormValue(currentForm);
    double minVal = minMaxForm.first;
    double maxVal = minMaxForm.second;
    formColor.clear();

    if (currentFormName == DUAL_0_FORM) {
        // Assign colors per dual vertex
        formColor.resize(mesh_dual->nVertices());
        for (Face f : mesh->faces()) {
            formColor[f.getIndex()] = mapToColor(currentForm[f.getIndex()], minVal, maxVal, "coolwarm");
        }
        for (Halfedge he : mesh->exteriorHalfedges()) {
            // extra vertices inserted in the middle of boundary edges
            formColor[heMap[he]] = mapToColor(currentForm[he.twin().face().getIndex()], minVal, maxVal, "coolwarm");
            // extra vertex per primal boundary vertex
            double avgVal = 0.0;
            double n = 0;
            for (Face f : he.vertex().adjacentFaces()) {
                avgVal += currentForm[f.getIndex()];
                n += 1.0;
            }
            avgVal /= n;
            formColor[vMap[he.vertex()]] = mapToColor(avgVal, minVal, maxVal, "coolwarm");
        }
    } else if (currentFormName == DUAL_1_FORM) {
        // 1-forms visualized with vectors, so make background white (via vertex
        // color)
        for (Vertex v : mesh_dual->vertices()) {
            formColor.push_back({1.0, 1.0, 1.0});
        }
        // Round very small values to 0 to avoid funky visualizations
        if (maxVal == 1e-5) {
            for (Edge e : mesh->edges()) {
                currentForm[e.getIndex()] = 0;
            }
        }
    } else if (currentFormName == DUAL_2_FORM) {
        // Assign colors per dual face
        formColor.resize(mesh_dual->nFaces());
        for (Vertex v : mesh->vertices()) {
            double A = geometry->barycentricDualArea(v);
            std::array<double, 3> color = mapToColor(currentForm[v.getIndex()] / A, minVal, maxVal, "coolwarm");
            // there is a 1-1 correspondence between primal vertices and dual
            // faces
            formColor[v.getIndex()] = color;
        }
    }
}

/*
 * Check if mesh is planar.
 */
bool isPlanar() {
    for (Vertex v : mesh->vertices()) {
        if (abs(geometry->inputVertexPositions[v].z) > 1e-8) {
            return false;
        }
    }
    return true;
}

/*
 * Visualize the current form.
 */
void showForm() {

    if (currentFormName == PRIMAL_0_FORM || currentFormName == PRIMAL_1_FORM || currentFormName == PRIMAL_2_FORM) {
        updateColorsPrimalForm();
        psMesh->setEnabled(true);
        psMesh_dual->setEnabled(false);
    } else {
        updateColorsDualForm();
        psMesh->setEnabled(false);
        psMesh_dual->setEnabled(true);
    }

    // color per primal vertex
    if (currentFormName == PRIMAL_0_FORM || currentFormName == PRIMAL_1_FORM) {
        polyscope::SurfaceVertexColorQuantity* formColorVertex = psMesh->addVertexColorQuantity("form", formColor);
        formColorVertex->setEnabled(true);
    }
    // color per dual vertex
    else if (currentFormName == DUAL_0_FORM || currentFormName == DUAL_1_FORM) {
        polyscope::SurfaceVertexColorQuantity* formColorVertex = psMesh_dual->addVertexColorQuantity("form", formColor);
        formColorVertex->setEnabled(true);
    }
    // color per primal face
    else if (currentFormName == PRIMAL_2_FORM) {
        polyscope::SurfaceFaceColorQuantity* formColorFace = psMesh->addFaceColorQuantity("form", formColor);
        formColorFace->setEnabled(true);
    }
    // color per dual face
    else if (currentFormName == DUAL_2_FORM) {
        polyscope::SurfaceFaceColorQuantity* formColorFace = psMesh_dual->addFaceColorQuantity("form", formColor);
        formColorFace->setEnabled(true);
    }

    // Visualize 1-form as a vector field if needed
    if (currentFormName == PRIMAL_1_FORM) {
        formColorEdge = psMesh->addFaceVectorQuantity("1form", interpolateWhitney());
        formColorEdge->setEnabled(true);
        formColorEdge->setVectorRadius(0.002);
        formColorEdge->setVectorLengthScale(lengthScalePrimal * 0.5);
        formColorEdge->setVectorColor({0.0, 0.0, 0.0});
    } else if (currentFormName == DUAL_1_FORM) {
        formColorDualEdge = psMesh_dual->addFaceVectorQuantity("1form", interpolateWachspressWhitney());
        formColorDualEdge->setEnabled(true);
        formColorDualEdge->setVectorRadius(0.002);
        formColorDualEdge->setVectorLengthScale(lengthScaleDual * 1.2);
        formColorDualEdge->setVectorColor({0.0, 0.0, 0.0});
    }
}


void redraw() {
    showForm(); // update form colors
    polyscope::requestRedraw();
}

/*
 * Generate random form.
 */
void randomize() {
    std::srand(std::time(0));
    if (currentFormName == PRIMAL_0_FORM || currentFormName == DUAL_2_FORM) {
        currentForm = generateRandomFormOnVertices(currentFormName == DUAL_2_FORM);
    } else if (currentFormName == PRIMAL_2_FORM || currentFormName == DUAL_0_FORM) {
        currentForm = generateRandomFormOnFaces(currentFormName == PRIMAL_2_FORM);
    } else {
        currentForm = generateRandomFormOnEdges();
    }
}

/*
 * User-defined buttons
 */
void functionCallback() {

    if (ImGui::Button("Randomize")) {
        randomize();
        redraw();
    }

    if (ImGui::Button("d")) {
        // Determine the operator and the next form.
        SparseMatrix<double> DEC_OP;
        int nextFormName;
        switch (currentFormName) {
        case PRIMAL_0_FORM:
            DEC_OP = D0;
            nextFormName = PRIMAL_1_FORM;
            break;
        case PRIMAL_1_FORM:
            DEC_OP = D1;
            nextFormName = PRIMAL_2_FORM;
            formColorEdge->setEnabled(false);
            break;
        case PRIMAL_2_FORM:
            // do nothing
            DEC_OP.resize(mesh->nFaces(), mesh->nFaces());
            DEC_OP.setIdentity();
            nextFormName = PRIMAL_2_FORM;
            break;
        case DUAL_0_FORM:
            DEC_OP = D1T;
            nextFormName = DUAL_1_FORM;
            break;
        case DUAL_1_FORM:
            DEC_OP = D0T;
            nextFormName = DUAL_2_FORM;
            formColorDualEdge->setEnabled(false);
            break;
        case DUAL_2_FORM:
            // do nothing
            DEC_OP.resize(mesh->nVertices(), mesh->nVertices());
            DEC_OP.setIdentity();
            nextFormName = DUAL_2_FORM;
            break;
        }
        currentForm = DEC_OP * currentForm;
        currentFormName = nextFormName;
        redraw();
    }

    if (ImGui::Button("*")) {
        // Determine the operator and the next form.
        SparseMatrix<double> DEC_OP;
        int nextFormName;
        switch (currentFormName) {
        case PRIMAL_0_FORM:
            DEC_OP = H0;
            nextFormName = DUAL_2_FORM;
            break;
        case PRIMAL_1_FORM:
            DEC_OP = H1;
            nextFormName = DUAL_1_FORM;
            formColorEdge->setEnabled(false);
            break;
        case PRIMAL_2_FORM:
            DEC_OP = H2;
            nextFormName = DUAL_0_FORM;
            break;
        case DUAL_0_FORM:
            DEC_OP = H2inv;
            nextFormName = PRIMAL_2_FORM;
            break;
        case DUAL_1_FORM:
            DEC_OP = H1inv;
            nextFormName = PRIMAL_1_FORM;
            formColorDualEdge->setEnabled(false);
            break;
        case DUAL_2_FORM:
            DEC_OP = H0inv;
            nextFormName = PRIMAL_0_FORM;
            break;
        }
        currentForm = DEC_OP * currentForm;
        currentFormName = nextFormName;
        redraw();
    }

    ImGui::Text("Current form: %s", formNameMap[currentFormName].c_str());

    showDiagram();
}


int main(int argc, char** argv) {
    // Configure the argument parser
    args::ArgumentParser parser("15-458 HW1");
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
    std::string filepath = "../../../input/hexagon.obj";
    if (inputFilename) {
        filepath = args::get(inputFilename);
    }

    // Load mesh
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(filepath);

    if (!isPlanar()) {
        std::cerr << "Please load a planar mesh with zero Z component" << std::endl;
        return EXIT_FAILURE;
    }

    // Initialize polyscope
    polyscope::view::style = polyscope::view::NavigateStyle::Planar;
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = functionCallback;

    // Add primal mesh to GUI
    psMesh = polyscope::registerSurfaceMesh("Primal mesh", geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));
    psMesh->setEnabled(true);
    lengthScalePrimal = geometry->meanEdgeLength();

    // Add visualization
    flipZ();
    createDualMesh();
    lengthScaleDual = geometry_dual->meanEdgeLength();
    // Add dualmesh to GUI
    psMesh_dual = polyscope::registerSurfaceMesh("Dual mesh", geometry_dual->inputVertexPositions,
                                                 mesh_dual->getFaceVertexList(), polyscopePermutations(*mesh_dual));
    psMesh_dual->setEnabled(false);

    // Build DEC operators.
    D0 = geometry->buildExteriorDerivative0Form();
    D1 = geometry->buildExteriorDerivative1Form();
    D0T = D0.transpose();
    D1T = D1.transpose();
    H0 = geometry->buildHodgeStar0Form();
    H1 = geometry->buildHodgeStar1Form();
    H2 = geometry->buildHodgeStar2Form();
    H0inv = sparseInverseDiagonal(H0);
    H1inv = sparseInverseDiagonal(H1);
    H2inv = sparseInverseDiagonal(H2);
    randomize();
    showForm();

    // initialize to something
    formColorEdge = psMesh->addFaceVectorQuantity("1form", std::vector<Vector3>(mesh->nFaces()));
    formColorDualEdge = psMesh_dual->addFaceVectorQuantity("1form", std::vector<Vector3>(mesh_dual->nFaces()));

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}