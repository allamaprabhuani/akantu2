Python examples
---------------

To run simulations in Python using Akantu, you first need to import the module:

.. code-block:: python

   import akantu as aka

The functions in Python are mostly the same as in C++.

The initiation differs. While in C++ a function is initialized using ``initialize("material.dat", argc, argv);``, in Python you should do:

.. code-block:: python

   aka.parseInput('material.dat')

The creation and loading of the mesh is done with:

.. code-block:: python

    mesh = aka.Mesh(spatial_dimension)
    mesh.read('mesh.msh')


.. include:: examples/python/solid_mechanics_model/README.rst
.. include:: examples/python/solid_mechanics_cohesive_model/README.rst
.. include:: examples/python/contact_mechanics_model/README.rst
.. include:: examples/python/phase_field_model/README.rst
