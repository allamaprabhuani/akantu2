Diffusion Model (2D and 3D)
```````````````````````````

:Sources:

   .. collapse:: heat_diffusion_static_2d.cc (click to expand)

      .. literalinclude:: examples/c++/diffusion_model/heat_diffusion_static_2d.cc
         :language: c++
         :lines: 20-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/c++/diffusion_model/material.dat
         :language: text

:Location:

   ``examples/c++/`` `diffusion_model <https://gitlab.com/akantu/akantu/-/blob/master/examples/c++/diffusion_model>`_


In ``diffusion_model``, examples of the ``HeatTransferModel`` are presented.

An example of a static heat propagation is presented in 
``heat_diffusion_static_2d.cc``. This example consists of a square 2D plate of 
:math:`1 \text{m}^2` having an initial temperature of :math:`100 \text{K}` 
everywhere but a non centered hot point maintained at 
:math:`300 \text{K}`. :numref:`fig-ex-diffusion_static` presents the geometry
of this case (left) and the results (right). The material used is a linear 
fictitious elastic material with a density of :math:`8940 \text{kg}/\text{m}^3`, 
a conductivity of :math:`401 \text{W}/\text{m}/\text{K}` and a specific heat 
capacity of :math:`385 \text{J}/\text{K}/\text{kg}`. 

The simulation is set following the procedure described in :ref:`sect-dm-using`

.. _fig-ex-diffusion_static:
.. figure:: examples/c++/diffusion_model/images/diffusion_static.png
            :align: center
            :width: 70%

            Initial (left) and final (right) temperature field 
            

:Sources:

   .. collapse:: heat_diffusion_dynamic_2d.cc (click to expand)

      .. literalinclude:: examples/c++/diffusion_model/heat_diffusion_static_2d.cc
         :language: c++
         :lines: 20-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/c++/diffusion_model/material.dat
         :language: text

In ``heat_diffusion_dynamics_2d.cc``, the same example is solved dynamically 
using an explicit time scheme. The time step used is :math:`0.12 \text{s}`. The only main difference with the previous example lies in the model initiation::

   model.initFull(_analysis_method = _explicit_lumped_mass);

.. _fig-ex-diffusion_explicit:
.. figure:: examples/c++/diffusion_model/images/hot-point-2.png
   :align: center     
   :width: 60%      
   
   Temperature field after 15000 time steps = 30 minutes. The lines represent 
   iso-surfaces.
   
:Sources:

   .. collapse:: heat_diffusion_dynamic_2d.cc (click to expand)

      .. literalinclude:: examples/c++/diffusion_model/heat_diffusion_static_2d.cc
         :language: c++
         :lines: 20-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/c++/diffusion_model/material.dat
         :language: text

In ``heat_diffusion_dynamics_3d.cc``, a 3D explicit dynamic heat propagation
problem is solved. It consists of a cube having an initial temperature of
:math:`100 \text{K}` everywhere but a centered sphere maintained at 
:math:`300 \text{K}`. 
The simulation is set exactly as ``heat_diffusion_dynamics_2d.cc`` except that the mesh is now a 3D mesh and that the heat source has a third coordinate and is placed at the cube center.
The mesh is initialized with::
   
   Int spatial_dimension = 3;
   Mesh mesh(spatial_dimension);
   mesh.read("cube.msh");

:numref:`fig-ex-diffusion_3d` presents the resulting temperature field evolution.
   
.. _fig-ex-diffusion_3d:
.. figure:: examples/c++/diffusion_model/images/diffusion_3d.gif
   :align: center     
   :width: 70%      
   
   Temperature field evolution.
