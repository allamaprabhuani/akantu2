solver_callback (2D)
''''''''''''''''''''

:Sources:

   .. collapse:: solver_callback.py (click to expand)

      .. literalinclude:: examples/python/solid_mechanics_model/solver_callback/solver_callback.py
         :language: python
         :lines: 9-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/python/solid_mechanics_model/solver_callback/material.dat
         :language: text

:Location:

   ``examples/python/solid_mechanics_model/`` `solver_callback <https://gitlab.com/akantu/akantu/-/blob/master/examples/python/solid_mechanics_model/solver_callback/>`_


In ``solver_callback``, it is shown how to write a solver callback function. It is done by writing a class that inherit 
from ``InterceptSolverCallback``::

	class SolverCallback(aka.InterceptSolverCallback):
		[.....]
	
	solver_callback = SolverCallback(model)
	
	model.solveStep(solver_callback)
	
In that case, the solver callback is used to modify the stiffness matrix.


