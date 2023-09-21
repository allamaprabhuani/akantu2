.. include:: auto_examples/index.rst


Tutorials
=========

Examples
=========

In addition to the tutorials, other C++ and Python examples are available in the ``examples`` folder. 

C++ examples
---------------------

Solid Mechanics Model
`````````````````````

Solid mechanics model examples are available in the ``solid_mechanics_model`` folder.

boundary_conditions
'''''''''''''''''''

In ``predifined_bc`` it is shown how to impose Dirichlet boundary condition using the predefined ``BC::Dirichlet::FixedValue`` (:numref:`fig-ex-predefined_bc`). Three built-in Dirichlet functors exist: ``FixedValue``, ``FlagOnly`` and ``IncrementValue``.

.. _fig-ex-predefined_bc:
.. figure:: figures/examples/predefined_bc.svg
            :align: center

            Dirichlet boundary conditions for the predifined_bc case. 

To define another functor, a class inherited from ``BC::Dirichlet::DirichletFunctor`` can be created as illustrated in the example ``user_defined_bc`` where a sinusoidal BC is imposed. The corresponding sinusoidal displacement is depicted in :numref:`fig-ex-user_defined_bc`. Note that a Neumann BC is also imposed.

.. _fig-ex-user_defined_bc:
.. figure:: figures/examples/user_defined_bc_displ_mag.png
            :align: center

            Displacement magnitude for the user_defined_bc example.
            
static
''''''

In ``static``, an example of how to solve a static problem is presented. The problem geometry is shown in :numref:`fig-ex-static`. The left and bottom side of a 2D plate are blocked with rollers and nodes from the left side are displaced upward by :math:`0.01\%`
of the length of the plate.

.. _fig-ex-static:
.. figure:: figures/examples/static_BC.svg
            :align: center

            Boundary conditions for the static example.
            
The solution for the static analysis is shown in :numref:`fig-ex-static_disp`.

.. _fig-ex-static_disp:
.. figure:: figures/examples/static_dipl_mag.svg
            :align: center

            Solution of the static analysis: displacement magnitude.
