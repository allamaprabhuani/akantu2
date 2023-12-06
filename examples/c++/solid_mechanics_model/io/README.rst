parser
''''''

:Sources:

   .. collapse:: example_parser.cc (click to expand)

      .. literalinclude:: examples/c++/solid_mechanics_model/io/parser/example_parser.cc
         :language: c++
         :lines: 20-

   .. collapse:: input_file.dat (click to expand)

      .. literalinclude:: examples/c++/solid_mechanics_model/io/parser/input_file.dat
         :language: text

:Location:

   ``examples/c++/solid_mechanics_model/io/`` `parser <https://gitlab.com/akantu/akantu/-/blob/master/examples/c++/solid_mechanics_model/io/parser>`_

In ``io/parser``, an example illustrating how to parse an input file with user-defined parameters is shown. As before, the text input file of the simulation is precised using the method ``initialize``. In the input file, additionally to the usual ``material elastic`` section, there is a section ``user parameters`` with extra user-defined parameters.
Within the main function, those parameters are retrived with::

   const ParserSection & usersect = getUserParser();
   Real parameter_name = usersect.getParameter("parameter_name");

dumper
''''''

:Sources:

   .. collapse:: dumper_low_level.cc (click to expand)

      .. literalinclude:: examples/c++/solid_mechanics_model/io/dumper/dumper_low_level.cc
         :language: c++
         :lines: 20-

   .. collapse:: dumpable_interface.cc (click to expand)

      .. literalinclude:: examples/c++/solid_mechanics_model/io/dumper/dumpable_interface.cc
         :language: c++
         :lines: 20-

:Location:

   ``examples/c++/solid_mechanics_model/io/`` `dumper <https://gitlab.com/akantu/akantu/-/blob/master/examples/c++/solid_mechanics_model/io/dumper>`_

In ``io/dumper``, examples of advanced dumping are shown.

``dumper_low_level`` aims at illustrating how to manipulate low-level methods of ``DumperIOHelper``. The goal is to visualize a colorized moving train with Paraview.
It is shown how to dump only a part of the mesh (here the wheels) using the function ``createElementGroup`` of the mesh object::

   ElementGroup & wheels_elements = mesh.createElementGroup("wheels", spatial_dimension);

One can then add an element to the group with::

   wheels_elements.append(mesh.getElementGroup("lwheel_1"));

where ``lwheel_1`` is the name of the element group in the mesh.

.. _fig-ex-train:
.. figure:: examples/c++/solid_mechanics_model/io/images/train.gif
            :align: center
            :width: 70%

            The wheels and the full train are dumped separately.

``dumpable_interface`` does the same as ``dumper_low_level`` but using ``dumpers::Dumpable`` which is an interface for other classes (Model, Mesh, ...) to dump themselves.
