parallel (2D)
'''''''''''''

:Sources:

   .. collapse:: parallel_2d.cc (click to expand)

      .. literalinclude:: examples/c++/solid_mechanics_model/parallel/parallel_2d.cc
         :language: c++
         :lines: 20-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/c++/solid_mechanics_model/parallel/material.dat
         :language: text

:Location:

   ``examples/c++/solid_mechanics_model/`` `parallel <https://gitlab.com/akantu/akantu/-/blob/master/examples/c++/solid_mechanics_model/parallel>`_

In ``parallel``, an example of how to run a 2D parallel simulation is presented.

First, the mesh needs to be distributed among the processors. This is done with::
    
    const auto & comm = Communicator::getStaticCommunicator();
    Int prank = comm.whoAmI();
    if (prank == 0) {
        mesh.read("square_2d.msh");
    }
    mesh.distribute();

Here ``prank`` designate the processor rank. The mesh is read only by processor 0 and then distributed among all the processors. 
All the prints are done only by processor 0. For instance::

    if (prank == 0) {
        std::cout << model.getMaterial(0) << std::endl;
    } 

.. figure:: examples/c++/solid_mechanics_model/parallel/images/parallel.png
            :align: center
            :width: 60%

            Displacement in the x direction.


