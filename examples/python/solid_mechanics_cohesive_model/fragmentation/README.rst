fragmentation (2D)
''''''''''''''''''

``fragmentation`` shows an example of a 1D bar fragmentation with extrinsic cohesive elements. It uses a custom boundary
condition to impose a constant velocity. This is done by creating a class that inherits from ``DirichletFunctor``. 
The result is shown in :numref:`fig-ex-fragmentation`. 

.. _fig-ex-fragmentation:
.. figure:: examples/python/solid_mechanics_cohesive_model/fragmentation/images/fragmentation.gif
            :align: center
            :width: 100%

            1D bar fragmentation.
            

