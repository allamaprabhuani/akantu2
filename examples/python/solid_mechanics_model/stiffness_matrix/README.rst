stiffness_matrix (2D)
'''''''''''''''''''''

:Sources:

   .. collapse:: stiffness_matrix.py (click to expand)

      .. literalinclude:: examples/python/solid_mechanics_model/stiffness_matrix/stiffness_matrix.py
         :language: python
         :lines: 9-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/python/solid_mechanics_model/stiffness_matrix/material.dat
         :language: text

:Location:

   ``examples/python/solid_mechanics_model/`` `stiffness_matrix <https://gitlab.com/akantu/akantu/-/blob/master/examples/python/solid_mechanics_model/stiffness_matrix/>`_


``stiffness_matrix`` shows how to get the stiffness matrix from a mesh. It is done with::
    model.assembleStiffnessMatrix()
    K = model.getDOFManager().getMatrix('K')
    stiff = aka.AkantuSparseMatrix(K).toarray()


