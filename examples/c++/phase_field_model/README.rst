Phase-field model
`````````````````

static
''''''

`phase_field_notch.cc` shows hot to set a quasi-static fracture simulation with phase-field. The geometry of the solid is a square plate with a notch. The loading is an imposed displacement at the top of the plate (mode I). 

.. figure:: examples/python/phase_field_model/images/phasefield-static-geo.svg
            :align: center
            :width: 50%

            Notched plate with boundary conditions and imposed displacement.

In addition, this example shows how to extract clusters delimited by the mesh boundaries and elements with damage values beyond a threshold. This can be useful to extract fragments after crack propagation.

.. figure:: examples/c++/phase_field_model/images/phase_field_clusters.png
            :align: center
            :width: 90%

            Damage field after a few iterations and two clusters (fragments) extracted.
