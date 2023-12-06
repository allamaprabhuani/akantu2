cohesive_extrinsic_ig_tg (2D)
'''''''''''''''''''''''''''''

:Sources:

   .. collapse:: cohesive_extrinsic_ig_tg.cc (click to expand)

      .. literalinclude:: examples/c++/solid_mechanics_cohesive_model/cohesive_extrinsic_ig_tg/cohesive_extrinsic_ig_tg.cc
         :language: c++
         :lines: 20-

   .. collapse:: material.dat (click to expand)

      .. _mat-cohesive-extrinsic-ig-tg:
      .. literalinclude:: examples/c++/solid_mechanics_cohesive_model/cohesive_extrinsic_ig_tg/material.dat
         :language: text
         :caption:

:Location:

   ``examples/c++/solid_mechanics_cohesive_model/`` `cohesive_extrinsic_ig_tg <https://gitlab.com/akantu/akantu/-/blob/master/examples/c++/solid_mechanics_cohesive_model/cohesive_extrinsic_ig_tg/>`_


In ``cohesive_extrinsic_ig_tg``, the insertion of cohesive element is not
limited to a given location. Rather, elements at the boundaries of the block and
those on the inside have a different critical stress ``sigma_c``. This is done
by defining two different materials in the :ref:`mat-cohesive-extrinsic-ig-tg`.
In this case the cohesive materials are chosen based on the bulk element on both
side. This is achieved by defining ``MaterialCohesiveRules``

The four block sides are then moved outwards. The resulting displacement is
shown in :numref:`fig-ex-cohesive-ext-ig-tg`.

.. _fig-ex-cohesive-ext-ig-tg:
.. figure:: examples/c++/solid_mechanics_cohesive_model/cohesive_extrinsic_ig_tg/images/cohesive_extrinsic_ig_tg.gif
            :align: center
            :width: 60%

            Displacement magnitude for the cohesive_extrinsic_ig_tg example.
