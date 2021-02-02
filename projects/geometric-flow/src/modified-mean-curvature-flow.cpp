// Implement member functions for ModifiedMeanCurvatureFlow class.
#include "modified-mean-curvature-flow.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
ModifiedMeanCurvatureFlow::ModifiedMeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry, and A (Laplace matrix)
    mesh = inputMesh;
    geometry = inputGeo;

    // TODO: build the Laplace matrix
    this->A = identityMatrix<double>(1); // placeholder
}

/*
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
SparseMatrix<double> ModifiedMeanCurvatureFlow::buildFlowOperator(const SparseMatrix<double>& M, double h) const {
    // TODO
    return identityMatrix<double>(1); // placeholder
}