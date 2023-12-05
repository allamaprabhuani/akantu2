stiffness_matrix
''''''''''''''''

``stiffness_matrix`` shows how to get the stiffness matrix from a mesh. It is done with::
    model.assembleStiffnessMatrix()
    K = model.getDOFManager().getMatrix('K')
    stiff = aka.AkantuSparseMatrix(K).toarray()


