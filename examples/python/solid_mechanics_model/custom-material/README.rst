custom-material
'''''''''''''''

In ``custom-material.py`` it is shown how to create a custom material behaviour. In this example, a linear elastic 
material is recreated. It is done by creating a class that inherits from ``aka.Material`` and register it 
to ``MaterialFactory``::

    class LocalElastic(aka.Material):
        [...]

    def allocator(_dim, unused, model, _id):
        return LocalElastic(model, _id)

    mat_factory = aka.MaterialFactory.getInstance()
    mat_factory.registerAllocator("local_elastic", allocator)
    
Wave propagation of a pulse in a bar fixed on the top, bottom and right boundaries is simulated using an explicit 
scheme. Results are shown in :numref:`fig-ex-plate_bar_custom`.

.. _fig-ex-plate_bar_custom:
.. figure:: examples/python/solid_mechanics_model/custom-material/images/pulse_bar_custom.gif
            :align: center
            :width: 90%

            Wave propagation in a bar.
            
In ``bi-material.py``, the same principle is used to create a bimaterial square. The displacement is shown in :numref:`fig-ex-square_custom`.

.. _fig-ex-square_custom:
.. figure:: examples/python/solid_mechanics_model/custom-material/images/square_displ.png
            :align: center
            :width: 60%

            Bimaterial square.
            

