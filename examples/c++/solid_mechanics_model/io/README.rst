parser
''''''

In ``io/parser``, an example illustrating how to parse an input file with user-defined parameters is shown. As before, the text input file of the simulation is precised using the method ``initialize``. In the input file, additionally to the usual ``material elastic`` section, there is a section ``user parameters`` with extra user-definied parameters.
Within the main function, those parameters are retrived with::

   const ParserSection & usersect = getUserParser();
   Real parameter_name = usersect.getParameter("parameter_name");

dumper
''''''

In ``io/dumper``, examples of advanced dumping are shown.

``dumpable_low_level`` aims at illustrating how to manipulate low-level methods of ``DumperIOHelper``. The goal is to visualize a colorized moving train with Paraview.
It is shown how to dump only a part of the mesh (here the wheels).

.. _fig-ex-train:
.. figure:: examples/c++/solid_mechanics_model/io/images/train.gif
            :align: center
            :width: 90%

            The wheels and the full train are dumped separately.

``dumpable_interface`` does the same as ``dumpable_low_level`` but using ``dumpers::Dumpable`` which is an interface for other classes (Model, Mesh, ...) to dump themselves.
