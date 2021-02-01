#include "distortion.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


/* Constructor
 * Input: The surface mesh <surfaceMesh> and geometry <geo>.
 */
Distortion::Distortion(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo) {
    mesh = surfaceMesh;
    geometry = geo;
}

/*
 * Compute the area-weighted quasi-conformal error of the given face.
 * See "geometric stretch metric" in http://hhoppe.com/ssp.pdf for more details.
 *
 * Input: A Face in the original mesh, and the parameterization.
 * Returns: The quasi conformal error of the face.
 */
double Distortion::computeQuasiConformalErrorPerFace(Face f, VertexData<Vector2>& flattening) {

    // Get vertices of face.
    Halfedge he = f.halfedge();
    Vertex x1 = he.vertex();
    Vertex x2 = he.next().vertex();
    Vertex x3 = he.next().next().vertex();

    // Edge vectors of original face
    Vector3 u1 = geometry->inputVertexPositions[x2] - geometry->inputVertexPositions[x1];
    Vector3 u2 = geometry->inputVertexPositions[x3] - geometry->inputVertexPositions[x1];
    // Edge vectors of face in parameterized mesh
    Vector2 v1 = flattening[x2] - flattening[x1];
    Vector2 v2 = flattening[x3] - flattening[x1];

    // Compute orthonormal bases.
    Vector3 e1 = u1.normalize();
    Vector3 e2 = u2 - (e1 * dot(u2, e1)); // subtract parallel part
    e2 = e2.normalize();

    Vector2 f1 = v1.normalize();
    Vector2 f2 = v2 - (f1 * dot(v2, f1));
    f2 = f2.normalize();

    // Compute vertex positions of each face projected onto basis.
    Vector2 p1 = {0.0, 0.0};
    Vector2 p2 = {dot(u1, e1), dot(u1, e2)};
    Vector2 p3 = {dot(u2, e1), dot(u2, e2)};

    Vector2 q1 = {0.0, 0.0};
    Vector2 q2 = {dot(v1, f1), dot(v1, f2)};
    Vector2 q3 = {dot(v2, f1), dot(v2, f2)};

    // Compute the entries of the Jacobian of the parameterization
    double A = 2.0 * cross(u1, u2).norm();
    Vector2 Ss = {0.0, 0.0};
    Ss += q1 * (p2[1] - p3[1]);
    Ss += q2 * (p3[1] - p1[1]);
    Ss += q3 * (p1[1] - p2[1]);
    Vector2 St = {0.0, 0.0};
    St += q1 * (p3[0] - p2[0]);
    St += q2 * (p1[0] - p3[0]);
    St += q3 * (p2[0] - p1[0]);
    Ss /= A;
    St /= A;

    // Compute the ratio of larger to smaller eigenvalues
    double a = dot(Ss, Ss);
    double b = dot(Ss, St);
    double c = dot(St, St);
    double det = sqrt((a - c) * (a - c) + 4.0 * b * b);
    double Gamma = sqrt(0.5 * (a + c + det));
    double gamma = sqrt(0.5 * (a + c - det));
    if (Gamma <= gamma) {
        return 1.0; // mapping is conformal
    }
    return Gamma / gamma; // ratio is > 1 otherwise
}

/*
 * Compute the average quasi-conformal error of the mesh resulting from a
 * parameterization algorithm.
 *
 * Input: The mesh parameterization, and a vector of colors per vertex (to be updated.)
 * Returns: The quasi conformal error of the mesh.
 */
double Distortion::computeQuasiConformalError(std::vector<std::array<double, 3>>& colors,
                                              VertexData<Vector2>& flattening) {

    double totalArea = geometry->totalArea(); // perhaps slightly more efficient to just compute in below loop
    double totalQcError = 0.0;
    std::vector<double> errorsPerFace(mesh->nFaces());

    for (Face f : mesh->faces()) {
        double qcError = computeQuasiConformalErrorPerFace(f, flattening);
        totalQcError += qcError * geometry->faceArea(f);
        errorsPerFace[f.getIndex()] = qcError;
    }

    // Assign colors per vertex for smoother visualization
    colors.resize(mesh->nVertices());
    for (Vertex v : mesh->vertices()) {
        double qcError = 0.0;
        for (Face f : v.adjacentFaces()) {
            qcError += errorsPerFace[f.getIndex()];
        }
        qcError /= v.faceDegree();
        colors[v.getIndex()] = hsv((2.0 - 4.0 * (qcError - 1.0)) / 3.0, 0.7, 0.65);
    }
    return totalQcError / totalArea;
}

/*
 * Compute the area scaling of the given face.
 *
 * Input: A Face in the original mesh, and the parameterization.
 * Returns: The area scaling of the face.
 */
double Distortion::computeAreaScalingPerFace(Face f, VertexData<Vector2>& flattening) {

    // Get vertices of face.
    Halfedge he = f.halfedge();
    Vertex p1 = he.vertex();
    Vertex p2 = he.next().vertex();
    Vertex p3 = he.next().next().vertex();

    // Edge vectors of original face
    Vector3 u1 = geometry->inputVertexPositions[p2] - geometry->inputVertexPositions[p1];
    Vector3 u2 = geometry->inputVertexPositions[p3] - geometry->inputVertexPositions[p1];
    // Edge vectors of face in parameterized mesh
    Vector2 v1 = flattening[p2] - flattening[p1];
    Vector2 v2 = flattening[p3] - flattening[p1];

    double u_area = cross(u1, u2).norm(); // 0.5 coefficients drop out
    double v_area = cross(v1, v2);

    return log(v_area / u_area);
}

/*
 * Compute the average area scaling of the mesh resulting from a
 * parameterization algorithm.
 *
 * Input: The mesh parameterization, and a vector of colors per vertex (to be updated.)
 * Returns: The total area scaling of the mesh.
 */
double Distortion::computeAreaScaling(std::vector<std::array<double, 3>>& colors, VertexData<Vector2>& flattening) {

    double totalArea = 0.0;
    double totalScaling = 0.0;
    FaceData<double> scalings = FaceData<double>(*mesh);
    double minScaling = std::numeric_limits<double>::infinity();
    double maxScaling = -std::numeric_limits<double>::infinity();
    double meanScaling = 0.0;
    std::vector<double> sortScalings(mesh->nFaces());

    for (Face f : mesh->faces()) {
        double scaling = computeAreaScalingPerFace(f, flattening);
        double area = geometry->faceArea(f);
        totalScaling += scaling * area;
        scalings[f] = scaling;
        meanScaling += scaling;
        sortScalings[f.getIndex()] = scaling;
        minScaling = std::min(minScaling, scaling);
        maxScaling = std::max(maxScaling, scaling);
    }
    if (minScaling == maxScaling) {
        minScaling = 0;
        maxScaling = 1e-3;
    }

    // Try to balance colors s.t. it is centered at 0.5
    // (seems to be mostly a problem on face.obj, where some faces near the corners of the eye throw the viz off)
    // Reset limits s.t. mean value coincides with mean value of colormap limits.
    int numFaces = mesh->nFaces();
    meanScaling /= numFaces;
    std::sort(sortScalings.begin(), sortScalings.end());
    // Heuristic to try to still keep as many values as possible: adjust the limit nearest to the tail of the
    // distribution.
    if (meanScaling < sortScalings[numFaces / 2]) {
        minScaling = 2 * meanScaling - maxScaling;
    } else {
        maxScaling = 2 * meanScaling - minScaling;
    }

    // Assign colors per vertex for smoother visualization
    colors.resize(mesh->nVertices());
    for (Vertex v : mesh->vertices()) {
        double scaling = 0.0;
        for (Face f : v.adjacentFaces()) {
            scaling += scalings[f];
        }
        scaling /= v.faceDegree();
        colors[v.getIndex()] = mapToColor(scaling, minScaling, maxScaling, "seismic");
    }

    return totalScaling / totalArea;
}