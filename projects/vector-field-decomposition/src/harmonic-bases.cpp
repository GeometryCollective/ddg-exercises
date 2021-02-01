// Implement member functions for HarmonicBases class.
#include "harmonic-bases.h"

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HarmonicBases::HarmonicBases(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build a closed, but not exact, primal 1-form ω.
 *
 * Input: A std::vector of Halfedges representing a homology generator of the mesh.
 * Returns: A vector representing a closed primal 1-form.
 */
Vector<double> HarmonicBases::buildClosedPrimalOneForm(const std::vector<Halfedge>& generator) const {

    // TODO
    return Vector<double>::Zero(1); // placeholder
}

/*
 * Compute the harmonic bases [γ1, γ2 ... γn] of the input mesh.
 *
 * Input: A std::vector of homology generators of the mesh (which are in turn represented as std::vectors of halfedges),
 * and a HodgeDecomposition object. Returns:
 */
std::vector<Vector<double>> HarmonicBases::compute(const std::vector<std::vector<Halfedge>>& generators,
                                                   const HodgeDecomposition& hodgeDecomposition) const {

    // TODO
    std::vector<Vector<double>> gammas;
    return gammas; // placeholder
}