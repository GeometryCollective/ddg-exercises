# ddg-exercises

This repo contains C++ skeleton code for course assignments from [Discrete Differential Geometry](https://brickisland.net/DDGSpring2020/) (15-458/858). 

For the JavaScript version, see https://github.com/cmu-geometry/ddg-exercises-js.

This code framework uses [Geometry Central](https://github.com/nmwsharp/geometry-central) for geometry processing utilities and [Polyscope](https://github.com/nmwsharp/polyscope) for visualization, which were developed by Nick Sharp and others in the [Geometry Collective](http://geometry.cs.cmu.edu/). Extensive documentation for these libraries ---_and how to build them on various platforms_--- can be found at the preceding links.  If you're having trouble building, please make sure to take a look before bugging the TAs! :-)  (We are of course still very happy to help if you're still having trouble.)

Documentation for Geometry Central can be found [here](https://geometry-central.net/).

Documentation for Polyscope can be found here [here](https://polyscope.run/).

## Getting started

Clone the repository and its submodules.
```
git clone --recursive https://github.com/GeometryCollective/ddg-exercises
cd ddg-exercises/projects
```

Each project in `ddg-exercises/projects` builds its own executable when compiled. To 
run a particular project `<project>`, go to the `projects/<project>` directory. The basic process for compiling is as follows. First, make a `build` directory and compile using
```
mkdir build
cd build
cmake ..
make
```
This builds an executable `main` which can then be run using
```
bin/main <optional_path_to_a_mesh>
```

(See [Geometry Central: Building](https://geometry-central.net/build/building/) for additional compiler flag options.

## Dependencies (all included)

1. Geometry processing and linear algebra - [Geometry Central](https://github.com/nmwsharp/geometry-central), which in turn has dependencies on [Eigen](https://eigen.tuxfamily.org) and [Suitesparse](https://people.engr.tamu.edu/davis/suitesparse.html).

2. Visualization - [Polyscope](https://github.com/nmwsharp/polyscope)

3. Unit tests - [Google Test](https://github.com/google/googletest)


## Author

Nicole Feng

Email: nfeng@cs.cmu.edu

Rohan Sawhney (original JavaScript version)

Email: rohansawhney@cs.cmu.edu

This code is directly based off Rohan's original JavaScript framework, [ddg-exercises-js](https://github.com/cmu-geometry/ddg-exercises-js).

## License

[MIT](https://opensource.org/licenses/MIT)
