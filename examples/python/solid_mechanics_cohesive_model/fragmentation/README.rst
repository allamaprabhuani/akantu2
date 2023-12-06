fragmentation (2D)
''''''''''''''''''

:Sources:

   .. collapse:: fragmentation.py (click to expand)

      .. literalinclude:: examples/python/solid_mechanics_cohesive_model/fragmentation/fragmentation.py
         :language: python
         :lines: 9-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/python/solid_mechanics_cohesive_model/fragmentation/material.dat
         :language: text

:Location:

   ``examples/python/solid_mechanics_cohesive_model/`` `fragmentation <https://gitlab.com/akantu/akantu/-/blob/master/examples/python/solid_mechanics_cohesive_model/fragmentation/>`_


``fragmentation`` shows an example of a 1D bar fragmentation with extrinsic cohesive elements. It uses a custom boundary
condition to impose a constant velocity. This is done by creating a class that inherits from ``DirichletFunctor``. 
The result is shown in :numref:`fig-ex-fragmentation`. 

.. _fig-ex-fragmentation:
.. figure:: examples/python/solid_mechanics_cohesive_model/fragmentation/images/fragmentation.gif
            :align: center
            :width: 100%

            1D bar fragmentation.
            

