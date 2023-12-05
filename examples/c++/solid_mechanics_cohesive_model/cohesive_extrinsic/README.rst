cohesive_extrinsic (2D)
'''''''''''''''''''''''

In ``cohesive_extrinsic``, an example of extrinsic cohesive elements is shown. 
An extrinsic simulation is initialized by setting ``_is_extrinsic`` argument of ``model.initFull`` to ``true``::
    
    model.initFull(_analysis_method = _explicit_lumped_mass, _is_extrinsic = true);

Cohesive elements are inserted during the simulation but at a location restricted with::
    
    CohesiveElementInserter & inserter = model.getElementInserter();
    inserter.setLimit(_y, 0.30, 0.20);
    model.updateAutomaticInsertion();

During the simulation, stress has to be checked in order to insert cohesive elements where the stress criterion is reached. This check is performed by calling the method ``checkCohesiveStress`` before each step resolution::
        
    model.checkCohesiveStress();
    model.solveStep();

For this specific example, a displacement is imposed based on the elements location. The corresponding displacements is shown in :numref:`fig-ex-cohesive-ext`.

.. _fig-ex-cohesive-ext:
.. figure:: examples/c++/solid_mechanics_cohesive_model/cohesive_extrinsic/images/cohesive_extrinsic.gif
            :align: center
            :width: 90%

            Displacement in the y direction for the cohesive_extrinsic example.
