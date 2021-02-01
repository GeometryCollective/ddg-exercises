This directory contains skeleton code for "Assignment 6: Vector Field Decomposition and Design" of the Discrete Differential Geometry course (15-458) at CMU. To get started, run

```
mkdir build
cd build
cmake ..
make
bin/main ../../../input/<mesh_file>
```

Initially, the program will just display the given mesh. (If a mesh is not provided, the default mesh is displayed.) 

For Assignment 6, you have the option of implementing either Hodge Decomposition or Trivial Connections. If you choose to implement Trivial Connections, you need to implement the following functions, assuming the input mesh is a topological sphere:
1. In `{CURRENT}/src/trivial-connections.cpp`:
	
	* constructor
	* `buildPeriodMatrix`
	* `computeCoExactComponent`
	* `computeHarmonicComponent`
	* `computeConnections`

In the GUI, click vertices to select them; shift-click to de-select. You can change the density value of an individual vertex by (re)-selecting it -- the currently active vertex will be red -- and entering in a value in the "Singularity" box. Hit the "Compute" button whenever you want to re-compute the solution.

Unit tests are built along with the rest of the program. To run the unit tests, run
```
bin/test-connect
```
This will output the results of several test functions in the terminal.

If there are bugs, please contact Nicole Feng (nfeng@cs.cmu.edu).