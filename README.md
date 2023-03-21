# `Akantu`: Swiss-Made Open-Source Finite-Element Library

`Akantu` means a little element in Kinyarwanda, a Bantu language. From now on it
is also an open- source object-oriented library which has the ambition to be
generic and efficient.

# Building `Akantu`

## Dependencies

In order to compile `Akantu`  any compiler supporting fully C++14 should work.
In addition some libraries are required:

 - CMake (>= 3.5.1)
 - Boost (preprocessor and Spirit)
 - zlib
 - blas/lapack

For the python interface:

 - Python (>=3 is recommended)
 - pybind11 (if not present the build system will try to download it)

To run parallel simulations:

 - MPI
 - Scotch

To use the static or implicit dynamic solvers at least one of the following libraries is needed:

 - MUMPS (since this is usually compiled in static you also need MUMPS dependencies)
 - PETSc

To compile the tests and examples:

 - Gmsh
 - google-test (if not present the build system will try to download it)

### On `.deb` based systems

``` sh
> sudo apt install cmake libboost-dev zlib1g-dev liblapack-dev libblas-dev gmsh
# For parallel
> sudo apt install mpi-default-dev libmumps-dev 
# For sequential
> sudo apt install libmumps-seq-dev 
```

## Configuring and compilation


`Akantu` is a [CMake](https://cmake.org/) project, so to configure it, you can follow the usual way:

``` sh
  > cd akantu
  > mkdir build
  > cd build
  > ccmake ..
  [ Set the options that you need ]
  > make
  > make install

```

## Using the python interface


You can install ``Akantu`` using pip, this will install a pre-compiled version:

``` sh
  > pip install akantu
```

You can then import the package in a python script as:

``` python
  import akantu
```

The python API is similar to the C++ one. If you
encounter any problem with the python interface, you are welcome to do a merge
request or post an issue on [GitLab](https://gitlab.com/akantu/akantu/-/issues).
  
# Tutorials with the python interface

To help getting started, multiple tutorials using the python interface are
available as notebooks with pre-installed version of `Akantu` on Binder. The
following tutorials are currently available:

[Plate whith a hole loaded](https://mybinder.org/v2/git/https%3A%2F%2Fgitlab.com%2Fakantu%2Ftutorials.git/HEAD?filepath=plate-hole/plate-hole.ipynb)

[Loaded cohesive crack](https://mybinder.org/v2/git/https%3A%2F%2Fgitlab.com%2Fakantu%2Ftutorials.git/HEAD?filepath=cohesive-fracture/cohesive-fracture.ipynb)

[Making your constitutive law in python](https://mybinder.org/v2/git/https%3A%2F%2Fgitlab.com%2Fakantu%2Ftutorials.git/HEAD?filepath=constitutive-laws/python_constitutive_law.ipynb)
