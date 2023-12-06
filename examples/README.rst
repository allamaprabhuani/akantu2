Examples
========

This part of the documentation describes the examples that are found in the
`examples <https://gitlab.com/akantu/akantu/-/tree/master/examples>`_ in the
repository. Akantu's example are separated in 2 types, the C++ and the python
ones, respectively in the 2 sub-folders ``c++/`` and ``python/``. The structure of
this documentation follows the folder structure of the ``examples`` folder.

The examples can be compiled by setting the option ``AKANTU_EXAMPLES`` in the
cmake configuration. They can be executed from after compilation from the
``examples`` folder in the build directory. Even though the python examples do not
need to be compiled, some file may still be generated it is then easier to run
them from the build directory. In order to set the different environment
variables needed a script ``akantu_environment.sh`` can be found in the build
directory.

Examples in both 2D and 3D are presented with the dimension is specified in the 
respective example titles. The only distinctions between a 2D and a 3D simulation lie in 
the mesh declaration. 
In C++::

    const Int spatial_dimension = 2;  // or 3 for 3D
    Mesh mesh(spatial_dimension);
    mesh.read("example_mesh.msh");

In Python::
    
    spatial_dimension = 2  # or 3 for 3D
    mesh = aka.Mesh(spatial_dimension)
    mesh.read("example_mesh.msh")

where ``example_mesh.msh`` is either a 2D or a 3D mesh.

.. include:: examples/c++/README.rst

.. include:: examples/python/README.rst
