Boundary conditions usage (2D)
''''''''''''''''''''''''''''''

In ``predifined_bc`` it is shown how to impose Dirichlet boundary condition
using the predefined ``BC::Dirichlet::FixedValue``
(:numref:`fig-ex-predefined_bc`). Three built-in Dirichlet functors exist:
``FixedValue``, ``FlagOnly`` and ``IncrementValue``.

.. _fig-ex-predefined_bc:
.. figure:: examples/c++/solid_mechanics_model/boundary_conditions/images/predefined_bc.svg
            :align: center

            Dirichlet boundary conditions for the predifined_bc case.

To define another functor, a class inherited from
``BC::Dirichlet::DirichletFunctor`` can be created as illustrated in the example
``user_defined_bc`` where a sinusoidal BC is imposed. The corresponding
sinusoidal displacement is depicted in :numref:`fig-ex-user_defined_bc`. Note
that a Neumann BC is also imposed.

.. _fig-ex-user_defined_bc:
.. figure:: examples/c++/solid_mechanics_model/boundary_conditions/images/user_defined_bc_displ_mag.png
            :align: center
            :width: 60%

            Displacement magnitude for the user_defined_bc example.
