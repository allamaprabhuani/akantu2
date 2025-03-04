cohesive_intrinsic (2D)
'''''''''''''''''''''''

:Sources:

   .. collapse:: cohesive_intrinsic.cc (click to expand)

      .. literalinclude:: examples/c++/solid_mechanics_cohesive_model/cohesive_intrinsic/cohesive_intrinsic.cc
         :language: c++
         :lines: 20-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/c++/solid_mechanics_cohesive_model/cohesive_intrinsic/material.dat
         :language: text

:Location:

   ``examples/c++/solid_mechanics_cohesive_model/`` `cohesive_intrinsic <https://gitlab.com/akantu/akantu/-/blob/master/examples/c++/solid_mechanics_cohesive_model/cohesive_intrinsic/>`_


In ``cohesive_intrinsic``, an example of intrinsic cohesive elements is shown. 
An intrinsic simulation is initialized by setting the ``_is_extrinsic`` argument of ``model.initFull`` to ``false``::
    
    model.initFull(_analysis_method = _explicit_lumped_mass, _is_extrinsic = false);

The cohesive elements are inserted between :math:`x = -0.26` and :math:`x =
-0.24` before the start of the simulation with
``model.getElementInserter().setLimit(_x, -0.26, -0.24);``. Elements to the
right of this limit are moved to the right. The resulting displacement is shown
in :numref:`fig-smm-cohesive-intrinsic`.

With intrinsic cohesive elements, a bi-linear or exponential cohesive law should
be used instead of a linear one (see section
:ref:`sect-smm-intrinsic-insertion`). This is set in the file ``material.dat``::

    material cohesive_bilinear [
	 name = cohesive
	 sigma_c = 1
	 beta = 1.5
	 G_c = 1
	 delta_0 = 0.1
	 penalty = 1e10
    ]

.. _fig-smm-cohesive-intrinsic:
.. figure:: examples/c++/solid_mechanics_cohesive_model/cohesive_intrinsic/images/cohesive_intrinsic.png
            :align: center
            :width: 60%

            Displacement in the x direction for the cohesive_intrinsic example.
