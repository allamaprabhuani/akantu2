cohesive_extrinsic_ig_tg
''''''''''''''''''''''''

In ``cohesive_extrinsic_ig_tg``, the insertion of cohesive element is not limited to a given location. Rather, elements at the boundaries of the block and those on the inside have a different critical stress ``sigma_c``. This is done by defining two different materials in the ``material.dat`` file::

    material cohesive_linear [
   	  name = tg_cohesive
   	  beta = 0
   	  G_c = 10
      sigma_c = 100
   ]

   material cohesive_linear [
      name = ig_cohesive
      beta = 0
      G_c = 10
      sigma_c = 20
   ]    


The four block sides are then moved outwards. The resulting displacement is shown in :numref:`fig-ex-cohesive-ext-ig-tg`.

.. _fig-ex-cohesive-ext-ig-tg:
.. figure:: examples/c++/solid_mechanics_cohesive_model/cohesive_extrinsic_ig_tg/images/cohesive_extrinsic_ig_tg.gif
            :align: center
            :width: 60%

            Displacement magnitude for the cohesive_extrinsic_ig_tg example.
