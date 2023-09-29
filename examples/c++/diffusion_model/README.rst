diffusion_model
'''''''''''''''

In ``diffusion_model``, examples of the ``HeatTransferModel`` are presented.

An example of a static heat propagation is presented in 
``heat_diffusion_static_2d.cc``. This example consists of a square 2D plate of 
:math:`1 \text{m}^2` having an initial temperature of :math:`100 \text{K}` 
everywhere but a none centered hot point maintained at 
:math:`300 \text{K}`. :numref:`fig-ex-diffusion_static` presents the geometry
of this case (left) and the results (right). The material used is a linear 
fictitious elastic material with a density of :math:`8940 \text{kg}/\text{m}^3`, 
a conductivity of :math:`401 \text{W}/\text{m}/\text{K}` and a specific heat 
capacity of :math:`385 \text{J}/\text{K}/\text{kg}`. 

.. _fig-ex-diffusion_static:
.. figure:: examples/c++/diffusion_model/images/diffusion_static.png
            :align: center
            :width: 70%

            Initial (left) and final (right) temperature field 
            

In ``heat_diffusion_dynamics_2d.cc``, the same example is solved dynamically 
using an explicit time scheme. The time step used is :math:`0.12 \text{s}`.

.. _fig-ex-diffusion_explicit:
.. figure:: examples/c++/diffusion_model/images/hot-point-1.png
   :align: center     
   :width: 70%      
   
   Temperature field after 15000 time steps = 30 minutes. The lines represent 
   iso-surfaces.
   
In ``heat_diffusion_dynamics_3d.cc``, a 3D explicit dynamic heat propagation
problem is solved. It consists of a cube having an initial temperature of
:math:`100 \text{K}` everywhere but a centered sphere maintained at 
:math:`300 \text{K}`. :numref:`fig-ex-diffusion_3d` presents the resulting 
temperature field evolution.
   
  .. _fig-ex-diffusion_3d:
.. figure:: examples/c++/diffusion_model/images/diffusion_3d.gif
   :align: center     
   :width: 70%      
   
   Temperature field evolution.
