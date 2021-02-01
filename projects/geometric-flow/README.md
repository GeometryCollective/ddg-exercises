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

2. In `{CURRENT}/src/mean-curvature-flow.cpp`:

	* `buildFlowOperator`
	* `integrate`

3. In `{CURRENT}/src/modified-mean-curvature-flow.cpp`:

	* constructor
	* `buildFlowOperator`

Once these functions are implemented, you can hit the "Mean curvature flow" or the 
"Modified MCF" buttons to apply their respective curvature flow operators with the given timestep.

Unit tests are built along with the rest of the program. To run the unit tests, run
```
bin/test-flow
```
This will output the results of several test functions in the terminal.

If there are bugs, please contact Nicole Feng (nfeng@cs.cmu.edu).
