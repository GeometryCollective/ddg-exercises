This directory contains skeleton code for "Assignment 5: Geodesic Distance" of the Discrete Differential Geometry course (15-458) at CMU. To get started, run

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
1. In `{CURRENT}/src/heat-method.cpp`:

	* constructor
	* `computeVectorField`
	* `computeDivergence`
	* `compute`

In the GUI, click vertices to select them; shift-click to de-select. The 
click-and-select functionality is currently not restricted to vertices, so 
selecting vertices in particular is a little cumbersome (sorry about that.)
Once the above functions are implemented, you can hit the "Solve" button to compute and 
visualize the solution.

Unit tests are built along with the rest of the program. To run the unit tests, run
```
bin/test-heat
```
This will output the results of several test functions in the terminal.


If there are bugs, please contact Nicole Feng (nfeng@cs.cmu.edu).
