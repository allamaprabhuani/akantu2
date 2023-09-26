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
