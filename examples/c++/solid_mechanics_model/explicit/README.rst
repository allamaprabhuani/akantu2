explicit (3D)
'''''''''''''

:Sources:

   .. collapse:: explicit_dynamic.cc (click to expand)

      .. literalinclude:: examples/c++/solid_mechanics_model/explicit/explicit_dynamic.cc
         :language: c++
         :lines: 20-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/c++/solid_mechanics_model/explicit/material.dat
         :language: text

:Location:

   ``examples/c++/solid_mechanics_model/`` `explicit <https://gitlab.com/akantu/akantu/-/blob/master/examples/c++/solid_mechanics_model/explicit>`_


In ``explicit``, an example of a 3D dynamic solution with an explicit time integration is shown.
The explicit scheme is selected using the ``_explicit_lumped_mass`` constant::

   model.initFull(_analysis_method = _explicit_lumped_mass);

Note that it is also the default value, hence using ``model.initFull();`` is equivalent.

This example models the propagation of a wave in a 3D steel beam. The beam and
the applied displacement in the :math:`x` direction are shown in
:numref:`fig-smm-explicit-bc`.

.. _fig-smm-explicit-bc:
.. figure:: examples/c++/solid_mechanics_model/explicit/images/explicit.svg
            :align: center
            :width: 90%

            Numerical setup.

The length and height of the beam are :math:`L={10}\mathrm{m}` and :math:`h =
{1}\mathrm{m}`, respectively. The material is linear elastic, homogeneous and
isotropic (density: :math:`{7800}\mathrm{kg}/\mathrm{m}^3`, Young's
modulus: :math:`{210}\mathrm{GPa}` and Poisson's ratio: :math:`0.3`). The
imposed displacement follow a Gaussian function with a maximum amplitude of
:math:`A = {0.01}\mathrm{m}`. The potential, kinetic and total energies are
computed. The safety factor is equal to :math:`0.8`.

The dynamic solution is depicted in :numref:`fig-smm-explicit-disp`.

.. _fig-smm-explicit-disp:
.. figure:: examples/c++/solid_mechanics_model/explicit/images/bar_pulse.gif
            :align: center
            :width: 100%

            Dynamic solution: lateral displacement.
