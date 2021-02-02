// Implement member functions for HodgeDecomposition class.
#include "hodge-decomposition.h"

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HodgeDecomposition::HodgeDecomposition(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;

    // TODO: build DEC operators
    this->hodge1 = identityMatrix<double>(1); // placeholder
    this->hodge2 = identityMatrix<double>(1); // placeholder
    this->d0 = identityMatrix<double>(1);     // placeholder
    this->d1 = identityMatrix<double>(1);     // placeholder

    // TODO: Build operator inverses.
    // Hint: Use the sparseInverseDiagonal() in utils/src/solvers.cpp to invert sparse diagonal matrices.
    this->hodge1Inv = identityMatrix<double>(1); // placeholder
    this->hodge2Inv = identityMatrix<double>(1); // placeholder
    this->d0T = identityMatrix<double>(1);       // placeholder
    this->d1T = identityMatrix<double>(1);       // placeholder

    // TODO: Construct 0-form Laplace matrix.
    // Shift matrix by a small constant (1e-8) to make it positive definite.
    this->A = identityMatrix<double>(1); // placeholder

    // TODO: Construct 2-form matrix.
    this->B = identityMatrix<double>(1); // placeholder
}

/*
 * Compute the 0-form potential α by solving the system 𝛿dα = 𝛿ω.
 *
 * Input: A primal 1-form on the edges of the input mesh.
 * Returns: The exact component dα of ω.
 */
Vector<double> HodgeDecomposition::computeExactComponent(const Vector<double>& omega) const {

    // TODO
    return Vector<double>::Zero(1); // placeholder
}

/*
 * Compute the 2-form potential β by solving the system d𝛿β = dω.
 *
 * Input: A primal 1-form on the edges of the input mesh.
 * Returns: The coexact component 𝛿β of ω.
 */
Vector<double> HodgeDecomposition::computeCoExactComponent(const Vector<double>& omega) const {

    // TODO
    return Vector<double>::Zero(1); // placeholder
}

/*
 * Compute the harmonic component γ = ω - dα - 𝛿β of ω.
 *
 * Input: A primal 1-form <omega> on the edges of the input mesh, the exact component <dAlpha> of ω, and the coexact
 * component <deltaBeta> of ω.
 * Returns: The coexact component 𝛿β of ω.
 */
Vector<double> HodgeDecomposition::computeHarmonicComponent(const Vector<double>& omega, const Vector<double>& dAlpha,
                                                            const Vector<double>& deltaBeta) const {

    // TODO
    return Vector<double>::Zero(1); // placeholder
}