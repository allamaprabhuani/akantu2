Solid Mechanics Cohesive Model
``````````````````````````````
Solid mechanics cohesive model examples are shown in ``solid_mechanics_cohesive_model``. This new model is called in a very similar way as the solid mechanics model::

   SolidMechanicsModelCohesive model(mesh);

Cohesive elements can be inserted intrinsically (when the mesh is generated) or extrinsically (during the simulation).

cohesive_intrinsic
''''''''''''''''''

In ``cohesive_intrinsic``, an example of intrinsic cohesive elements is shown. The cohesive elements are inserted between :math:`x = -0.26` and :math:`x = -0.24` before the start of the simulation with ``model.getElementInserter().setLimit(_x, -0.26, -0.24);``. Elements to the right of this limit are moved to the right. The resulting displacement is shown in :numref:`fig-ex-cohesive-int`.

.. _fig-ex-cohesive-int:
.. figure:: examples/c++/solid_mechanics_cohesive_model/cohesive_intrinsic/images/cohesive_intrinsic.png
            :align: center
            :width: 70%

            Displacement in the x direction for the cohesive_intrinsic example.
            
cohesive_extrinsic
''''''''''''''''''

In ``cohesive_extrinsic``, cohesive elements are inserted during the simulation but at a location restricted with::
    CohesiveElementInserter & inserter = model.getElementInserter();
    inserter.setLimit(_y, 0.30, 0.20);
    model.updateAutomaticInsertion();
A displacement is then imposed based on the elements location. The corresponding displacements is shown in :numref:`fig-ex-cohesive-ext`.

.. _fig-ex-cohesive-ext:
.. figure:: examples/c++/solid_mechanics_cohesive_model/cohesive_extrinsic/images/cohesive_extrinsic.gif
            :align: center
            :width: 90%

            Displacement in the y direction for the cohesive_extrinsic example.

cohesive_extrinsic_ig_tg
''''''''''''''''''''''''

In ``cohesive_extrinsic_ig_tg``, the insertion of cohesive element is not limited to a given location. Rather, elements at the boundaries of the block and those on the inside have a different critical stress defined in the ``material.dat`` file. The four block sides are then moved outwards. The resulting displacement is shown in :numref:`fig-ex-cohesive-ext-ig-tg`.

.. _fig-ex-cohesive-ext-ig-tg:
.. figure:: examples/c++/solid_mechanics_cohesive_model/cohesive_extrinsic_ig_tg/images/cohesive_extrinsic_ig_tg.gif
            :align: center
            :width: 100%

            Displacement magnitude for the cohesive_extrinsic_ig_tg example.
