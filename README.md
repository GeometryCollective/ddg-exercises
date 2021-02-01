# ddg-exercises-cpp

ddg-exercises-cpp is the C++ version of [ddg-exercises-js](https://github.com/cmu-geometry/ddg-exercises-js), which contains skeleton code for various geometry processing algorithms. In particular, it uses [Geometry Central](https://github.com/nmwsharp/geometry-central) for geometry processing utilities and [Polyscope](https://github.com/nmwsharp/polyscope) for visualization, which were developed by Nick Sharp and others in the [Geometry Collective](http://geometry.cs.cmu.edu/).

Documentation for Geometry Central can be found [here](https://geometry-central.net/).

Documentation for Polyscope can be found here [here](https://polyscope.run/).

## Getting started

Clone the repository and its submodules.
```
git clone --recursive <url-of-this-repo>
cd ddg-exercises-cpp/projects
```

Each project in `ddg-exercises-cpp/projects` builds its own executable when compiled. To 
run a particular project `<project>`, go to the `projects/<project>` directory. First, compile using
```
mkdir build
cd build
cmake ..
make
```
which builds an executable `main` which can then be run using
```
bin/main <optional_path_to_a_mesh>
```

## Dependencies (all included)

1. Geometry processing and linear algebra - [Geometry Central](https://github.com/nmwsharp/geometry-central), which in turn has dependencies on [Eigen](https://eigen.tuxfamily.org) and [Suitesparse](https://people.engr.tamu.edu/davis/suitesparse.html).

2. Visualization - [Polyscope](https://github.com/nmwsharp/polyscope)

3. Unit tests - [Google Test](https://github.com/google/googletest)


## Author

Nicole Feng

Email: nfeng@cs.cmu.edu

## License

[MIT](https://opensource.org/licenses/MIT)
