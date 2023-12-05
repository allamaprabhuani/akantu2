plate_hole (2D)
'''''''''''''''

In ``plate_hole``, an example of a static solution of a loaded plate with a hole. Because of the symmetries, only a quarter of the problem is modeled. The corresponding domain and boundary conditions is shown in :numref:`fig-ex-plate_hole`.
The static solver is initialized with::

    model.initFull(_analysis_method=aka._static)

Boundary conditions are applied with::

    model.applyBC(aka.FixedValue(0.0, aka._x), "XBlocked")
    model.applyBC(aka.FixedValue(0.0, aka._y), "YBlocked")
    model.applyBC(aka.FromTraction(trac), "Traction")

where ``trac`` is a numpy array of size ``spatial_dimension``.

Additionally, the simulation can be run in parallel. To do so, the following lines are added::
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        prank = comm.Get_rank()
    except ImportError:
        prank = 0

Similarly to C++, the mesh has to be distributed between the processors with::
    
    mesh.distribute()

.. _fig-ex-plate_hole:
.. figure:: examples/python/solid_mechanics_model/plate-hole/images/plate_hole.svg
            :align: center
            :width: 30%

            Plate with a hole geometry.
            
The displacement magnitude is displayed in :numref:`fig-ex-plate_hole_displ`.

.. _fig-ex-plate_hole_displ:
.. figure:: examples/python/solid_mechanics_model/plate-hole/images/plate_hole_displ_mag.png
            :align: center
            :width: 40%

            Displacement magnitude.

