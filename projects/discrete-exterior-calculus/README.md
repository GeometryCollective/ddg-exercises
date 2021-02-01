This directory contains skeleton code for "Assignment 1: Exterior Calculus" of the Discrete Differential Geometry course (15-458) at CMU. To get started, run

```
mkdir build
cd build
cmake ..
make
bin/main ../../../input/<mesh_file>
```

Initially, the program will just display the given mesh. (If a mesh is not 
provided, the default mesh is displayed.) 
You can view the wireframe for the primal mesh by checking the "Edges" box under 
the structure labeled "Primal mesh" on the left menu, and likewise for the dual mesh.

For this assignment, you need to implement the following functions:
1. In `{HOME}/core/src/geometry.cpp`:

	* `cotan`
	* `barycentricDualArea`

2. In `{HOME}/core/src/discrete-exterior-calculus.cpp`:

	* `buildHodgeStar0Form`
	* `buildHodgeStar1Form`
	* `buildHodgeStar2Form`
	* `buildExteriorDerivative0Form`
	* `buildExteriorDerivative1Form`

Unit tests are built along with the rest of the program. To run the unit tests, run
```
bin/test-dec
```
This will output the results of several test functions in the terminal.

If there are bugs, please contact Nicole Feng (nfeng@cs.cmu.edu).
