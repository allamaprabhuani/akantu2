Scalability test (3D)
`````````````````````

:Sources:

   .. collapse:: cohesive_extrinsic.cc (click to expand)

      .. literalinclude:: examples/c++/scalability_test/cohesive_extrinsic.cc
         :language: c++
         :lines: 20-

   .. collapse:: material-elastic.dat (click to expand)

      .. literalinclude:: examples/c++/scalability_test/material-elastic.dat
         :language: text

:Location:

   ``examples/c++/`` `scalability_test <https://gitlab.com/akantu/akantu/-/blob/master/examples/c++/scalability_test>`_

  This example is used to do scalability test, with elastic material or cohesive elements inserted on the fly.
  The `cube.geo` should generate a mesh with roughly 4'500'000 elements and 730'000 nodes. The `cube.msh` file included is a `tiny` cube only to test if the code works.

  To run the full simulation the full mesh has to be generated `gmsh -3 cube.geo -o cube.msh`

  The simulation consist of a cube with a compressive force on top and a shear force on four sides as shown in the figure bellow

  .. figure:: examples/c++/scalability_test/images/cube.svg
            :align: center
            :width: 60%

  This examples was used to run the scalability test from the JOSS paper the
  full mesh and results presented in the following graph can be found here:
  `perf-test-akantu-cohesive
  <https://gitlab.com/akantu/performance-testing-cohesive/-/tree/publications/joss?ref_type=tags>`_

  .. figure:: examples/c++/scalability_test/images/TTS.svg
            :align: center
            :width: 60%
