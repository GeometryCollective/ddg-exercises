This directory contains skeleton code for "Assignment 2: Investigating Curvature" of the Discrete Differential Geometry course (15-458) at CMU. To get started, run

```
mkdir build
cd build
cmake ..
make
bin/main ../../../input/<mesh_file>
```

Initially, the program will just display the given mesh. (If a mesh is not 
provided, the default mesh is displayed.) You will notice that you can toggle which quantities are displayed on the mesh via the menu on the left, but for convenience you can simply use the buttons provided in the callback UI in the upper right.

For this assignment, you need to implement the following functions:
1. In `{HOME}/core/src/geometry.cpp`:

	* `angle`
	* `dihedralAngle`
	* `vertexNormalAngleWeighted`
	* `vertexNormalSphereInscribed`
	* `vertexNormalAreaWeighted`
	* `vertexNormalGaussianCurvature`
	* `vertexNormalMeanCurvature`
	* `angleDefect`
	* `totalAngleDefect`
	* `scalarMeanCurvature`
	* `circumcentricDualArea`
	* `principalCurvatures`

Unit tests are built along with the rest of the program. To run the unit tests, run
```
bin/test-curv
```
This will output the results of several test functions in the terminal.

If there are bugs, please contact Nicole Feng (nfeng@cs.cmu.edu).
