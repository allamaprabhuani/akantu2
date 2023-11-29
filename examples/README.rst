Examples
========

This part of the documentation describes the examples that are found in the
`examples <https://gitlab.com/akantu/akantu/-/tree/master/examples>`_ in the
repository. Akantu's example are separated in 2 types, the C++ and the python
ones, respectively in the 2 sub-folders `c++/` and `python/`. The structure of
this documentation follows the folder structure of the `examples` folder.

The examples can be compiled by setting the option `AKANTU_EXAMPLES` in the
cmake configuration. They can be executed from after compilation from the
`examples` folder in the build directory. Even though the python examples do not
need to be compiled, some file may still be generated it is then easier to run
them from the build directory. In order to set the different environment
variables needed a script `akantu_environment.sh` can be found in the build
directory.


.. include:: examples/c++/README.rst

.. include:: examples/python/README.rst
