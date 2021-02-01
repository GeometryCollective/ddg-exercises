// Implement member functions for TrivialConnections class.
#include "trivial-connections.h"

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
TrivialConnections::TrivialConnections(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;

    // TODO: Build harmonic bases
    this->bases; // placeholder;

    // Build period matrix.
    this->P = this->buildPeriodMatrix();

    // TODO: Store DEC operators
    this->A = identityMatrix<double>(1);      // placeholder
    this->hodge1 = identityMatrix<double>(1); // placeholder
    this->d0 = identityMatrix<double>(1);     // placeholder
}

/*
 * Builds the period matrix Pij = ∑_{ek ∈ li} (ξj)k, where li is the ith homology generator, ek is a dual edge in li and
 * ξj is the jth harmonic 1-form basis.
 *
 * Input:
 * Returns: A sparse matrix represending the period matrix.
 */
SparseMatrix<double> TrivialConnections::buildPeriodMatrix() const {
    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Determine if a mesh satisfies Gauss-Bonnet.
 *
 * Input: A vector where the ith entry is the the index of the singularity at the ith vertex.
 * Returns: True if mesh satisfies Gauss-Bonnet, false otherwise.
 */
bool TrivialConnections::satsifyGaussBonnet(const Vector<double>& singularity) const {

    return (abs(singularity.sum() - geometry->eulerCharacteristic()) < 1e-8);
}

/*
 * Compute the dual 0-form potential β by solving the system d𝛿β = -K + 2π * singularity.
 *
 * Input: A vector where the ith entry is the the index of the singularity at the ith vertex.
 * Returns: The coexact component 𝛿β.
 */
Vector<double> TrivialConnections::computeCoExactComponent(const Vector<double>& singularity) const {

    // TODO
    return Vector<double>::Zero(1); // placeholder
}


/*
 * Given an initial angle αi in face i, this function computes the new angle αj in the neighboring face j as
 * αj = αi - θij + θji, where θij and θji are the angles between the shared edge e and an arbitrary but fixed reference
 * direction in faces i and j. Repeating this procedure for n consecutive dual edges in a generator gives a sequence of
 * angles α0, . . . , αn with a resulting total angle defect equal to αn - α0. This corresponds to transporting a vector
 * around a generator by unfolding, sliding and refolding it across neighboring faces without any extra in plane
 * rotation.
 *
 * Input: A halfedge lying on the shared edge between face i and j, and the initial angle αi.
 * Returns: The new angle αj.
 */
double TrivialConnections::transportNoRotation(Halfedge he, double alphaI) const {

    Vector3 u = geometry->halfedgeVector(he);

    // Compute two orthonormal tangent vectors for each face.
    Face fi = he.face();
    Face fj = he.twin().face();
    Vector3 e1 = geometry->halfedgeVector(fi.halfedge()).normalize();
    Vector3 e2 = cross(geometry->faceNormal(fi), e1);
    Vector3 f1 = geometry->halfedgeVector(fj.halfedge()).normalize();
    Vector3 f2 = cross(geometry->faceNormal(fj), f1);
    double thetaIJ = atan2(dot(u, e2), dot(u, e1));
    double thetaJI = atan2(dot(u, f2), dot(u, f1));

    return alphaI - thetaIJ + thetaJI;
}

/*
 * Compute the harmonic component γ = ∑_{i = 1, ..., 2g} zi ξi by solving the system Pz = v - ∑𝛿β.
 * v - ∑𝛿β should be normalized to lie between -π and π.
 *
 * Input: The coexact component 𝛿β.
 * Returns: The harmonic component γ.
 */
Vector<double> TrivialConnections::computeHarmonicComponent(const Vector<double>& deltaBeta) const {

    // TODO
    return Vector<double>::Zero(1); // placeholder
}

/*
 * Compute the dual 1-form connections φ = 𝛿β + γ.
 *
 * Input: A vector where the ith entry is the the index of the singularity at the ith vertex.
 * Returns: A vector representing the connections.
 */
Vector<double> TrivialConnections::computeConnections(const Vector<double>& singularity) const {

    if (!this->satsifyGaussBonnet(singularity)) {
        std::cerr << "Singularities do not add up to the Euler characteristic of the mesh" << std::endl;
        return Vector<double>::Zero(mesh->nEdges());
    }
    // TODO: Compute connections on topological spheres
    return Vector<double>::Zero(1); // placeholder
}