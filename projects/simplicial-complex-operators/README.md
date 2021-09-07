This directory contains skeleton code for "Assignment 0: Combinatorial Surfaces" of the Discrete Differential Geometry course (15-458) at CMU. To get started, run

```
mkdir build
cd build
cmake ..
make -j4
bin/main ../../../input/<mesh_file>
```

Initially, the program will just display the given mesh. (If a mesh is not  provided, the default mesh is displayed.) 

Turn on wireframe visualization by checking the "Edges" box on the lefthand side menu. Use the mouse to click and select 
vertices, edges, and faces on screen. The click-and-select functionality works best if you click near the center of each mesh element.

For this assignment, you need to implement the following functions:
1. In `{CURRENT}/src/simplicial-complex-operators.cpp`:

	* `assignElementIndices`
	* `buildVertexEdgeAdjacencyMatrix`
	* `buildEdgeFaceAdjacencyMatrix`
	* `buildVertexVector`
	* `buildEdgeVector`
	* `buildFaceVector`
	* `star`
	* `closure`
	* `link`
	* `isComplex`
	* `isPureComplex`
	* `boundary`

Once the functions are filled in correctly, the buttons in the GUI will apply the corresponding simplicial operators.
You can also adjust visualization settings of the mesh and other structures in the lefthand-side menus.

In addition, the file `core/include/mesh_subset.h` contains the implementation of the `MeshSubset` class, which has functions for constructing, duplicating, etc. subsets.

Unit tests are built along with the rest of the program. To run the unit tests, run
```
bin/test-sco
```
This will output the results of several test functions in the terminal.

If there are bugs, or if you simply have questions, please contact Nicole Feng (nfeng@cs.cmu.edu).
