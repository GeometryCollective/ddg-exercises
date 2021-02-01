This directory contains skeleton code for "Assignment 3: The Laplacian" of the Discrete Differential Geometry course (15-458) at CMU. To get started, run

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
1. In `{HOME}/core/src/geometry.cpp`:

	* `laplaceMatrix`
	* `massMatrix`

2. In `{CURRENT}/src/scalar-poisson-problem.cpp`:

	* constructor
	* `solve`

In the GUI, click vertices to select them; shift-click to de-select. 
The  click-and-select functionality is currently not restricted to vertices, so 
selecting vertices in particular is a little cumbersome (sorry about that.)
You can change the density value of an individual vertex by (re)-selecting it -- 
the currently active vertex will be red -- and entering in a value in the 
"Density" box. Hit the "Solve" button whenever you want to re-compute the solution.

If there are bugs, please contact Nicole Feng(nfeng @cs.cmu.edu).
