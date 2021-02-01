This directory contains skeleton code for "Assignment 4: Conformal Parameterization" of the Discrete Differential Geometry course (15-458) at CMU. To get started, run

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

	* `complexLaplaceMatrix`

2. In `{CURRENT}/src/spectral-conformal-parameterization.cpp`:

	* `buildConformalEnergy`
	* `flatten`

3. In `{CURRENT}/src/solvers.cpp`:

	* `solveInversePowerMethod`
	* `residual`

Once these functions are implemented, you will be able to toggle the original/flattened mesh in the GUI, as well as view different error functions on the mesh. You can also view the mesh wireframe by 
toggling the "Edges" option under the mesh structure in the lefthand-side menu.

Unit tests are built along with the rest of the program. To run the unit tests, run
```
bin/test-param
```
This will output the results of several test functions in the terminal.

If there are bugs, please contact Nicole Feng (nfeng@cs.cmu.edu).
