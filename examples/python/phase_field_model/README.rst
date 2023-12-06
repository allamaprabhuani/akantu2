Phase-field model
`````````````````

static
''''''

:Sources:

   .. collapse:: phasefield-static.py (click to expand)

      .. literalinclude:: examples/python/phase_field_model/phasefield-static.py
         :language: python
         :lines: 9-

   .. collapse:: material_static.dat (click to expand)

      .. literalinclude:: examples/python/phase_field_model/material_static.dat
         :language: text

:Location:

   ``examples/python/`` `phase_field_model <https://gitlab.com/akantu/akantu/-/blob/master/examples/python/phase_field_model>`_


`phasefield-static.py` shows how to setup a static phase-field fracture simulation. An imposed displacement is imposed on the top of a notched square plate.

.. figure:: examples/python/phase_field_model/images/phasefield-static-geo.svg
            :align: center
            :width: 50%

            Notched plate with boundary conditions and imposed displacement.

In static simulations, we use loading steps to apply the displacement incrementally. At each time step, the solvers of the solid mechanics model and the phase-field model are called alternately until convergence is reached.

.. figure:: examples/python/phase_field_model/images/phasefield-static.png
            :align: center
            :width: 50%

            Damage field after a few iterations.

dynamic
'''''''

:Sources:

   .. collapse:: phasefield-dynamic.py (click to expand)

      .. literalinclude:: examples/python/phase_field_model/phasefield-dynamic.py
         :language: python
         :lines: 9-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/python/phase_field_model/material.dat
         :language: text

`phasefield-dynamic.py` shows how to setup a dynamic phase-field fracture simulation. A notched plate is pre-strained in mode I using Dirichlet BC and a static solve. The simulation is then continued in dynamic using an explicit Neumark scheme.

.. figure:: examples/python/phase_field_model/images/phasefield-dynamic-geo.svg
            :align: center
            :width: 80%

            Notched plate with boundary conditions and imposed displacement.

At each time step, each solver is called once to find the displacement field and the damage field.

.. figure:: examples/python/phase_field_model/images/phasefield-dynamic.png
            :align: center
            :width: 80%

            Crack propagation and branching.
