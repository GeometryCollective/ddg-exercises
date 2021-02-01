This directory contains skeleton code for "Assignment 6: Vector Field Decomposition and Design" of the Discrete Differential Geometry course (15-458) at CMU. To get started, run

```
mkdir build
cd build
cmake ..
make
bin/main ../../../input/<mesh_file>
```

Initially, the program will just display the given mesh. (If a mesh is not 
provided, the default mesh is displayed.) 

For this project, you need to implement the following functions:
1. In `{CURRENT}/src/tree-cotree.cpp`:

	* `buildPrimalSpanningTree`
	* `inPrimalSpanningTree`
	* `buildDualSpanningCoTree`
	* `inDualSpanningCotree`
	* `buildGenerators`

2. In `{CURRENT}/src/harmonic-bases.cpp`:

	* `buildClosedPrimalOneForm`
	* `compute`

For Assignment 6, you have the option of implementing either Hodge Decomposition or Trivial Connections. If you choose to implement Hodge Decomposition, you need to implement the following functions:
1. In `{CURRENT}/src/hodge-decomposition.cpp`:
	* constructor
	* `computeExactComponent`
	* `computeCoExactComponent`
	* `computeHarmonicComponent`

You will notice that you can toggle which quantities are displayed on the mesh via the menu on the left, but for convenience you can simply use the buttons provided in the callback UI in the upper right. You can also adjust visualization settings of the mesh and other structures in the lefthand-side menus.

Unit tests are built along with the rest of the program. To run the unit tests, run
```
bin/test-decomp
```
This will output the results of several test functions in the terminal.

If there are bugs, please contact Nicole Feng (nfeng@cs.cmu.edu).
