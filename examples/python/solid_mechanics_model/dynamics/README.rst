dynamics (2D)
'''''''''''''

:Sources:

   .. collapse:: dynamics.py (click to expand)

      .. literalinclude:: examples/python/solid_mechanics_model/dynamics/dynamics.py
         :language: python
         :lines: 9-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/python/solid_mechanics_model/dynamics/material.dat
         :language: text

:Location:

   ``examples/python/solid_mechanics_model/`` `dynamics <https://gitlab.com/akantu/akantu/-/blob/master/examples/python/solid_mechanics_model/dynamics/>`_

In ``dynamics``, an example of a dynamic simulation solved with an explicit time intergration is shown. This examples 
model the propagation of a wave in a bar. The geometry is depicted in :numref:`fig-ex-dynamics_geom`. A pulse is impose 
as an initial displacement over the bar. Results are depicted in :numref:`fig-ex-dynamics_displ`.

The explicit time integration scheme is set with:

.. code-block:: python

    model.initFull(_analysis_method=aka._explicit_lumped_mass)

This example shows how to create a new functor to set the boundary conditions. This is done by creating a class that inherits from ``aka.FixedValue`` (``MyFixedValue(aka.FixedValue)`` in this case).
The boundary conditions are then applied with:

.. code-block:: python
    
    model.applyBC(MyFixedValue(0, aka._x), "XBlocked")
    model.applyBC(MyFixedValue(0, aka._y), "YBlocked")   

.. _fig-ex-dynamics_geom:
.. figure:: examples/python/solid_mechanics_model/dynamics/images/bar_geom.svg
            :align: center
            :width: 70%

            Plate with a hole geometry.

.. _fig-ex-dynamics_displ:
.. figure:: examples/python/solid_mechanics_model/dynamics/images/bar.gif
            :align: center
            :width: 70%

            Displacement magnitude.

