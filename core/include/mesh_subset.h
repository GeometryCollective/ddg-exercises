#pragma once

#include <set>

class MeshSubset {

  public:
    std::set<size_t> vertices;
    std::set<size_t> edges;
    std::set<size_t> faces;

    /* Initialize an empty MeshSubset. */
    MeshSubset() {}

    /* Initialize a MeshSubset with the given vertices, edges, and faces. */
    MeshSubset(const std::set<size_t>& V, const std::set<size_t>& E, const std::set<size_t>& F) {

        vertices = V;
        edges = E;
        faces = F;
    }

    /* Make a deep copy of the input MeshSubset and return it as a new
     * MeshSubset.
     */
    MeshSubset deepCopy() const {

        std::set<size_t> newVertices = vertices;
        std::set<size_t> newEdges = edges;
        std::set<size_t> newFaces = faces;
        return MeshSubset(newVertices, newEdges, newFaces);
    }

    /* Add a vertex to this subset. */
    void addVertex(size_t index) {
        vertices.insert(index);
    }

    /* Add a set of vertices to this subset. */
    void addVertices(const std::set<size_t>& V) {
        for (std::set<size_t>::iterator it = V.begin(); it != V.end(); ++it) {
            vertices.insert(*it);
        }
    }

    /* Delete a vertex from this subset. */
    void deleteVertex(size_t index) {
        vertices.erase(index);
    }

    /* Delete a set of vertices from this subset. */
    void deleteVertices(const std::set<size_t>& V) {
        for (std::set<size_t>::iterator it = V.begin(); it != V.end(); ++it) {
            vertices.erase(*it);
        }
    }

    /* Add an edge to this subset. */
    void addEdge(size_t index) {
        edges.insert(index);
    }

    /* Add a set of edges to this subset. */
    void addEdges(const std::set<size_t>& E) {
        for (std::set<size_t>::iterator it = E.begin(); it != E.end(); ++it) {
            edges.insert(*it);
        }
    }

    /* Delete an edge from this subset. */
    void deleteEdge(size_t index) {
        edges.erase(index);
    }

    /* Delete a set of edges from this subset. */
    void deleteEdges(const std::set<size_t>& E) {
        for (std::set<size_t>::iterator it = E.begin(); it != E.end(); ++it) {
            edges.erase(*it);
        }
    }

    /* Add a face to this subset. */
    void addFace(size_t index) {
        faces.insert(index);
    }

    /* Add a set of faces to this subset. */
    void addFaces(const std::set<size_t>& F) {
        for (std::set<size_t>::iterator it = F.begin(); it != F.end(); ++it) {
            faces.insert(*it);
        }
    }

    /* Delete a face from this subset. */
    void deleteFace(size_t index) {
        faces.erase(index);
    }

    /* Delete a set of faces from this subset. */
    void deleteFaces(const std::set<size_t>& F) {
        for (std::set<size_t>::iterator it = F.begin(); it != F.end(); ++it) {
            faces.erase(*it);
        }
    }

    /* Returns true if subsets are equivalent. */
    bool equals(const MeshSubset& other) const {
        // == compares elements at each position; but std::set always orders elements upon insertion/initialization
        return (vertices == other.vertices) && (edges == other.edges) && (faces == other.faces);
    }

    /* Adds a subset's vertices, edges, and faces to this subset. */
    void addSubset(const MeshSubset& other) {
        this->addVertices(other.vertices);
        this->addEdges(other.edges);
        this->addFaces(other.faces);
    }

    /* Removes a subset's vertices, edges, and faces from this subset. */
    void deleteSubset(const MeshSubset& other) {
        this->deleteVertices(other.vertices);
        this->deleteEdges(other.edges);
        this->deleteFaces(other.faces);
    }

    /* Print vertices. */
    void printVertices() {
        std::cout << "Vertices: ";
        for (std::set<size_t>::iterator it = vertices.begin(); it != vertices.end(); ++it) {
            std::cout << *it << ", ";
        }
        std::cerr << std::endl;
    }

    /* Print edges. */
    void printEdges() {
        std::cout << "Edges: ";
        for (std::set<size_t>::iterator it = edges.begin(); it != edges.end(); ++it) {
            std::cout << *it << ", ";
        }
        std::cerr << std::endl;
    }

    /* Print Faces. */
    void printFaces() {
        std::cout << "Faces: ";
        for (std::set<size_t>::iterator it = faces.begin(); it != faces.end(); ++it) {
            std::cout << *it << ", ";
        }
        std::cerr << std::endl;
    }
};