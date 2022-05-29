// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"
#include <Eigen/Sparse>
#include <vector>
#include <iostream>

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.

    typedef Eigen::Triplet<size_t> T;

    std::vector<T> triplets; 
    triplets.reserve(mesh->nEdges() * 2);

    for (Edge edge: mesh->edges()) {
        triplets.push_back(T(edge.getIndex(), edge.firstVertex().getIndex(), 1 ));
        triplets.push_back(T(edge.getIndex(), edge.secondVertex().getIndex(), 1 ));
    }

    SparseMatrix<size_t> a0(mesh->nEdges(), mesh->nVertices());
    a0.setFromTriplets(triplets.begin(), triplets.end());
    return a0;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    typedef Eigen::Triplet<size_t> T;

    std::vector<T> triplets; 
    triplets.reserve(mesh->nFaces() * 3);

    for (Face face: mesh->faces()) {
        for(Edge edge: face.adjacentEdges()) {
            triplets.push_back(T(face.getIndex(), edge.getIndex(), 1 ));
        }
    }

    SparseMatrix<size_t> a1(mesh->nFaces(), mesh->nEdges());
    a1.setFromTriplets(triplets.begin(), triplets.end());
    return a1;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    Vector<size_t> encoding = Vector<size_t>::Zero(mesh->nVertices());
    for (auto i: subset.vertices) {
        encoding[i] = 1;
    }
    return encoding;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    Vector<size_t> encoding = Vector<size_t>::Zero(mesh->nEdges());
    for (auto i: subset.edges) {
        encoding[i] = 1;
    }
    return encoding;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    Vector<size_t> encoding = Vector<size_t>::Zero(mesh->nFaces());
    for (auto i: subset.faces) {
        encoding[i] = 1;
    }
    return encoding;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    auto starVertices = buildVertexVector(subset);   
    auto starEdges = buildEdgeVector(subset);   
    auto starFaces = buildFaceVector(subset);   

    auto vertexToEdge = buildVertexEdgeAdjacencyMatrix();
    starEdges += vertexToEdge * starVertices;

    auto edgeToFace = buildFaceEdgeAdjacencyMatrix();
    starFaces += edgeToFace * starEdges;

    MeshSubset star = MeshSubset();
    star.addVertices(subset.vertices);

    for(Edge edge: mesh->edges()) {
        if (starEdges(edge.getIndex()) > 0) {
            star.addEdge(edge.getIndex());
        }
    }

    for(Face face: mesh->faces()) {
        if (starFaces(face.getIndex()) > 0) {
            star.addFace(face.getIndex());
        }
    }

    return star;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    auto closureVertices = buildVertexVector(subset);   
    auto closureEdges = buildEdgeVector(subset);   
    auto closureFaces = buildFaceVector(subset);   

    auto edgeToFace = buildFaceEdgeAdjacencyMatrix();
    closureEdges += edgeToFace.transpose() * closureFaces;

    auto vertexToEdge = buildVertexEdgeAdjacencyMatrix();
    closureVertices += vertexToEdge.transpose() * closureEdges;
    
    MeshSubset closure = MeshSubset();
    closure.addFaces(subset.faces);

    for(Edge edge: mesh->edges()) {
        if (closureEdges(edge.getIndex()) > 0) {
            closure.addEdge(edge.getIndex());
        }
    }

    for(Vertex vertex: mesh->vertices()) {
        if (closureVertices(vertex.getIndex()) > 0) {
            closure.addVertex(vertex.getIndex());
        }
    }

    return closure;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    auto outer = closure(star(subset));
    auto inner = star(closure(subset));

    std::set<size_t> vertices;
    std::set_difference(
        outer.vertices.begin(), outer.vertices.end(), 
        inner.vertices.begin(), inner.vertices.end(), 
        std::inserter(vertices, vertices.begin()));

    std::set<size_t> edges;
    std::set_difference(
        outer.edges.begin(), outer.edges.end(), 
        inner.edges.begin(), inner.edges.end(), 
        std::inserter(edges, edges.begin()));
    
    std::set<size_t> faces;

    MeshSubset link = MeshSubset(vertices, edges, faces);
    return link;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    auto closureSubset = closure(subset);
    return closureSubset.equals(subset);

    // auto abstractFaces = buildFaceVector(subset);   
    // auto abstractEdges = buildEdgeVector(subset);  
    // auto abstractVertices = buildVertexVector(subset);   
     
    // auto edgeToFace = buildFaceEdgeAdjacencyMatrix();
    // auto edges = edgeToFace.transpose() * abstractFaces;

    // for(Edge edge: mesh->edges()) {
    //     size_t i = edge.getIndex();
    //     if (edges(i) > 0 && subset.edges.find(i) == subset.edges.end() ) {
    //         return false;
    //     }
    // }
    
    // auto vertexToEdge = buildVertexEdgeAdjacencyMatrix();
    // auto vertices = vertexToEdge.transpose() * abstractEdges;

    // for(Vertex vertex: mesh->vertices()) {
    //     size_t i = vertex.getIndex();
    //     if (vertices(i) > 0 && subset.vertices.find(i) == subset.vertices.end() ) {
    //         return false;
    //     }
    // }

    // return true;
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    if ( !isComplex(subset) ) return -1;

    if (subset.vertices.size() == 0) return -1;

    if (subset.edges.size() == 0) return 0;

    if (subset.faces.size() == 0) return 1;

    auto complexEdges = buildEdgeVector(subset);  
    auto vertexToEdge = buildVertexEdgeAdjacencyMatrix();
    auto vertices = vertexToEdge.transpose() * complexEdges;

    for(Vertex vertex: mesh->vertices()) {
        size_t i = vertex.getIndex();
        if (subset.vertices.find(i) != subset.vertices.end() && vertices(i) == 0) {
            return -1;
        }
    }

    if (subset.faces.size() == 0) return 1;

    auto complexFaces = buildFaceVector(subset);   
    auto edgeToFace = buildFaceEdgeAdjacencyMatrix();
    auto edges = edgeToFace.transpose() * complexFaces;

    for(Edge edge: mesh->edges()) {
        size_t i = edge.getIndex();
        if (subset.edges.find(i) != subset.edges.end() && edges(i) == 0) {
            return -1;
        }
    }

    return 2;
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    auto complexFaces = buildFaceVector(subset);   
    auto edgeToFace = buildFaceEdgeAdjacencyMatrix();
    auto edges = edgeToFace.transpose() * complexFaces;

    MeshSubset boundary = MeshSubset();

    for(Edge edge: mesh->edges()) {
        size_t i = edge.getIndex();
        if (subset.edges.find(i) != subset.edges.end() && edges(i) == 1) {
            boundary.addEdge(i);
        }
    }

    auto complexEdges = buildEdgeVector(subset);   
    auto vertexToEdge = buildVertexEdgeAdjacencyMatrix();
    auto vertices = vertexToEdge.transpose() * complexEdges;

    for(Vertex vertex: mesh->vertices()) {
        size_t i = vertex.getIndex();
        if (subset.vertices.find(i) != subset.vertices.end() && vertices(i) == 1) {
            boundary.addVertex(i);
        }
    }

    boundary = closure(boundary);
    return boundary;
}