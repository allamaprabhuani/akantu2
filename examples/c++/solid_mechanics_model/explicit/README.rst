explicit
''''''''

Corresponding files:
 - `explicit_dynamic.cc <https://gitlab.com/akantu/akantu/-/blob/master/examples/c++/solid_mechanics_model/explicit/explicit_dynamic.cc>`_
 - `material.dat <https://gitlab.com/akantu/akantu/-/blob/master/examples/c++/solid_mechanics_model/explicit/material.dat>`_

In ``explicit``, an example of a dynamic solution with an explicit time integration is shown.
The explicit scheme is selected using the ``_explicit_lumped_mass`` constant::

   model.initFull(_analysis_method = _explicit_lumped_mass);

Note that it is also the default value, hence using ``model.initFull();`` is equivalent.

This example models the propagation of a wave in a steel beam. The beam and the
applied displacement in the :math:`x` direction are shown in
:numref:`fig-ex-explicit`.

.. _fig-ex-explicit:
.. figure:: examples/c++/solid_mechanics_model/explicit/images/explicit.svg
            :align: center
            :width: 90%

            Numerical setup.

The length and height of the beam are :math:`L=\SI{10}{\meter}` and :math:`h =
\SI{1}{\meter}`, respectively. The material is linear elastic, homogeneous and
isotropic (density: :math:`\SI{7800}{\kilo\gram\per\meter\cubed}`, Young's
modulus: :math:`\SI{210}{\giga\pascal}` and Poisson's ratio: :math:`0.3`). The
imposed displacement follow a Gaussian function with a maximum amplitude of
:math:`A = \SI{0.01}{\meter}`. The potential, kinetic and total energies are
computed. The safety factor is equal to :math:`0.8`.

The dynamic solution is depicted in :numref:`fig-ex-explicit_disp`.

.. _fig-ex-explicit_disp:
.. figure:: examples/c++/solid_mechanics_model/explicit/images/bar_pulse.gif
            :align: center
            :width: 100%

            Dynamic solution: lateral displacement.

.. literalinclude:: examples/c++/solid_mechanics_model/explicit/explicit_dynamic.cc
   :language: c++

