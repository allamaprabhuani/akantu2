new_material
''''''''''''

:Sources:

   .. collapse:: new_local_material.cc (click to expand)

      .. literalinclude:: examples/c++/solid_mechanics_model/new_material/new_local_material.cc
         :language: c++
         :lines: 20-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/c++/solid_mechanics_model/new_material/material.dat
         :language: text

:Location:

   ``examples/c++/solid_mechanics_model/`` `new_material   <https://gitlab.com/akantu/akantu/-/blob/master/examples/c++/solid_mechanics_model/new_material>`_


In ``new_material`` it is shown how to use a user-defined material for the simulation. All the details are described in :ref:`sect-smm-ncl`. The geometry solved is shown in :numref:`fig-ex-new_material`.

.. _fig-ex-new_material:
.. figure:: examples/c++/solid_mechanics_model/new_material/images/barre_trou.svg
            :align: center
            :width: 70%

            Problem geometry.
