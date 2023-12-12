Static test cases (2D)
''''''''''''''''''''''

:Sources:

   .. collapse:: static.cc (click to expand)

      .. literalinclude:: examples/c++/solid_mechanics_model/static/static.cc
         :language: c++
         :lines: 20-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/c++/solid_mechanics_model/static/material.dat
         :language: text

:Location:

   ``examples/c++/solid_mechanics_model/`` `static <https://gitlab.com/akantu/akantu/-/blob/master/examples/c++/solid_mechanics_model/static>`_


In ``static``, an example of how to solve a static problem is presented. The
problem geometry is shown in :numref:`fig-ex-static`. The left and bottom side
of a 2D plate are blocked with rollers and nodes from the left side are
displaced upward by :math:`0.01\%` of the length of the plate.

.. _fig-ex-static:
.. figure:: examples/c++/solid_mechanics_model/static/images/static_BC.svg
            :align: center

            Boundary conditions for the static example.

The solution for the static analysis is shown in :numref:`fig-ex-static_disp`.

.. _fig-ex-static_disp:
.. figure:: examples/c++/solid_mechanics_model/static/images/static_displ_mag.png
            :align: center
            :width: 60%

            Solution of the static analysis: displacement magnitude.
