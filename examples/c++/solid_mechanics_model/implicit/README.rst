implicit
''''''''

In ``implicit``, an example of a dynamic solution with an implicit time integration is shown.
The implicit scheme is selected using the ``_implicit_dynamic`` constant::

   model.initFull(_analysis_method = _implicit_dynamic);

This example consists of
a 3D beam of
:math:`10\mathrm{m}\times1\mathrm{m}\times1\mathrm{m}` blocked
on one side and is on a roller on the other side. A constant force of
:math:`5\mathrm{kN}` is applied in its middle.
:numref:`fig-ex-implicit-dynamic` presents the geometry of this case. The
material used is a fictitious linear elastic material with a density of
:math:`1000 \mathrm{kg/m}^3`, a Young's Modulus of
:math:`120 \mathrm{MPa}` and Poisson's ratio of :math:`0.3`. These values
were chosen to simplify the analytical solution.

An approximation of the dynamic response of the middle point of the
beam is given by:

.. math::

    u\left(\frac{L}{2}, t\right)
    = \frac{1}{\pi^4} \left(1 - cos\left(\pi^2 t\right) +
    \frac{1}{81}\left(1 - cos\left(3^2 \pi^2 t\right)\right) +
    \frac{1}{625}\left(1 - cos\left(5^2 \pi^2 t\right)\right)\right)

.. _fig-ex-implicit-dynamic:
.. figure:: examples/c++/solid_mechanics_model/implicit/images/implicit_dynamic.svg
            :align: center
            :width: 75%

            Numerical setup.

..
   \begin{figure}[!htb]
     \centering
     \includegraphics[scale=.6]{figures/implicit_dynamic}
     \caption{Numerical setup}
     \label{fig-smm-implicit:dynamic}
   \end{figure}

:numref:`fig-ex-implicit-dynamic_solution` presents the deformed
beam at 3 different times during the simulation: time steps 0, 1000 and
2000.

.. _fig-ex-implicit-dynamic_solution:
.. figure:: examples/c++/solid_mechanics_model/implicit/images/dynamic_analysis.png
            :align: center
            :width: 40%

            Deformed beam at three different times (displacement :math:`\times
            10`).
