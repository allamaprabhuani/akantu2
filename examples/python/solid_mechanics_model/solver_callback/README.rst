solver_callback (2D)
''''''''''''''''''''

In ``solver_callback``, it is shown how to write a solver callback function. It is done by writing a class that inherit 
from ``InterceptSolverCallback``::

	class SolverCallback(aka.InterceptSolverCallback):
		[.....]
	
	solver_callback = SolverCallback(model)
	
	model.solveStep(solver_callback)
	
In that case, the solver callback is used to modify the stiffness matrix.


