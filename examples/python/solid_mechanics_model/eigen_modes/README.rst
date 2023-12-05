eigen_modes (2D)
''''''''''''''''

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

            
