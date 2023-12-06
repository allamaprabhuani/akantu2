eigen_modes (2D)
''''''''''''''''

:Sources:

   .. collapse:: eigen_modes.py (click to expand)

      .. literalinclude:: examples/python/solid_mechanics_model/eigen_modes/eigen_modes.py
         :language: python
         :lines: 9-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/python/solid_mechanics_model/eigen_modes/material.dat
         :language: text

:Location:

   ``examples/python/solid_mechanics_model/`` `eigen_modes <https://gitlab.com/akantu/akantu/-/blob/master/examples/python/solid_mechanics_model/eigen_modes/>`_


In ``eigen_modes`` it is shown how to compute the eigen modes using the library `scipy.sparse`. The mode can be specify with the ``-m`` argument. Simulation is performed on a bar where a pulse is imposed. 
The ``-p`` argument will plot the following figures :numref:`fig-ex-eigen`:

.. _fig-ex-eigen:
.. figure:: examples/python/solid_mechanics_model/eigen_modes/images/eigen_modes.png
            :align: center
            :width: 80%

            Energy norms as a fonction of time (left), space-time diagram for diplacements (center) and space-time 
            diagram for velocities (right) with the default values.

An implicit time integration scheme is used and set with::
    model.initFull(aka._implicit_dynamic)

            
