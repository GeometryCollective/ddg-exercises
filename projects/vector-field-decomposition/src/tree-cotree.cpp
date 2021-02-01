// Implement member functions for TreeCotree class.
#include "tree-cotree.h"

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
TreeCotree::TreeCotree(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build a primal spanning tree on a mesh without boundary. More specifically, populate the member variable
 * <vertexParent>, which is a std::map that maps each vertex of the input mesh to its parent in the primal spanning
 * tree.
 *
 * Input:
 * Returns:
 */
void TreeCotree::buildPrimalSpanningTree() {

    // TODO
}

/*
 * Check whether a halfedge is in the primal spanning tree.
 *
 * Input: A halfedge <he>
 * Returns: True if <he> is in the primal spanning tree, false otherwise.
 */
bool TreeCotree::inPrimalSpanningTree(Halfedge he) {

    // TODO
    return false; // placeholder
}

/*
 * Build a dual spanning tree on a mesh without boundary. More specificially, populate the member variable <faceParent>,
 * which is a std::map that maps each face of the input mesh to its parent in the dual spanning tree.
 *
 * Input:
 * Returns:
 */
void TreeCotree::buildDualSpanningCoTree() {

    // TODO
}

/*
 * Check whether a halfedge is in the dual spanning tree.
 *
 * Input: A halfedge <he>
 * Returns: True if <he> is in the dual spanning tree, false otherwise.
 */
bool TreeCotree::inDualSpanningCotree(Halfedge he) {

    // TODO
    return false; // placeholder
}

/*
 * Returns a halfedge lying on the shared edge between face f and g.
 *
 * Input: Two adjacent faces <f> and <g>.
 * Returns: A halfedge lying on the shared edge between face f and g.
 */
Halfedge TreeCotree::sharedHalfedge(Face f, Face g) const {

    for (Halfedge he : f.adjacentHalfedges()) {
        if (he.twin().face() == g) {
            return he;
        }
    }
    // Should never get here!
    std::cerr << "Oops, TreeCotree::sharedHalfedge() received bad input." << std::endl;
    return f.halfedge();
}

/*
 * Compute the homology generators of the mesh.
 *
 * Input:
 * Returns:
 */
void TreeCotree::buildGenerators() {

    // order doesn't matter in a mesh without boundary
    buildPrimalSpanningTree();
    buildDualSpanningCoTree();

    // TODO: Build generators and populate this->generators
}