cohesive_intrinsic
''''''''''''''''''

In ``cohesive_intrinsic``, an example of intrinsic cohesive elements is shown. The cohesive elements are inserted between :math:`x = -0.26` and :math:`x = -0.24` before the start of the simulation with ``model.getElementInserter().setLimit(_x, -0.26, -0.24);``. Elements to the right of this limit are moved to the right. The resulting displacement is shown in :numref:`fig-ex-cohesive-int`.

.. _fig-ex-cohesive-int:
.. figure:: examples/c++/solid_mechanics_cohesive_model/cohesive_intrinsic/images/cohesive_intrinsic.png
            :align: center
            :width: 60%

            Displacement in the x direction for the cohesive_intrinsic example.
