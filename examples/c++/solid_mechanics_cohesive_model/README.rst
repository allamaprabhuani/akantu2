Solid Mechanics Cohesive Model
``````````````````````````````
Solid mechanics cohesive model examples are shown in ``solid_mechanics_cohesive_model``. This new model is called in a very similar way as the solid mechanics model::

   SolidMechanicsModelCohesive model(mesh);

Cohesive elements can be inserted intrinsically (when the mesh is generated) or extrinsically (during the simulation).

.. include:: examples/c++/solid_mechanics_cohesive_model/cohesive_intrinsic/README.rst

.. include:: examples/c++/solid_mechanics_cohesive_model/cohesive_extrinsic/README.rst

.. include:: examples/c++/solid_mechanics_cohesive_model/cohesive_extrinsic_ig_tg/README.rst
