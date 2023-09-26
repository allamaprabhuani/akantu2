cohesive
''''''''

In ``cohesive/plate.py``, an example of extrinsic cohesive elements insertion is shown. This example simulates a plate 
with a pre-existing crack pulled. The geometry is depicted in :numref:`fig-ex-cohesive_plate`. 

.. _fig-ex-cohesive_plate:
.. figure:: examples/python/solid_mechanics_cohesive_model/cohesive/images/plate.svg
            :align: center
            :width: 100%

            Problem geometry.
            
The problem is solved statically until the onset of instability and is then solved dynamically with an explicit
integration scheme. This is done by adding a new solver with ``model.initNewSolver()``::
    model = aka.SolidMechanicsModelCohesive(mesh)
    model.initFull(_analysis_method=aka._static, _is_extrinsic=True)
    model.initNewSolver(aka._explicit_lumped_mass)
    [....]
    model.solveStep("static")
    [....]
    model.solveStep("explicit_lumped")
    
The crack opening is displayed in :numref:`fig-ex-cohesive_plate_gif`.

.. _fig-ex-cohesive_plate_gif:
.. figure:: examples/python/solid_mechanics_cohesive_model/cohesive/images/plate.gif
            :align: center
            :width: 100%

            Stresses in the plate.

