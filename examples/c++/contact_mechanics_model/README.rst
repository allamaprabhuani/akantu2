Contact mechanics model (2D)
````````````````````````````

`contact_explicit_static` and `contact_explicit_dynamic` are solving a 2D Hertz contact patch test using the ``CouplerSolidContact``.
The two examples follow what is described extensively in section :ref:`_sect-cmm-coupling-with-smm`. The only main difference between `contact_explicit_static` and `contact_explicit_dynamic` is the solver used::
    
    // for contact_explicit_static
    coupler.initFull(_analysis_method = _static);  
    // for contact_explicit_dynamic
    coupler.initFull(_analysis_method = _explicit_lumped_mass);  

The ``material.dat`` file contain a ``contact_detector`` and a ``contact_resolution penalty_linear`` section as explains in section :ref:`_sect-cmm-contact-detection`.

.. figure:: examples/c++/contact_mechanics_model/images/hertz.svg
            :align: center
            :width: 60%

.. figure:: examples/c++/contact_mechanics_model/images/hertz.png
            :align: center
            :width: 60%
