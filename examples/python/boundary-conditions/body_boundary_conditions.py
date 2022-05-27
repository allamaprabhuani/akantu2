#!/usr/bin/env python3
""" plate.py: Plate with hole linear elastic example"""

__author__ = "Guillaume Anciaux"
__credits__ = [
    "Guillaume Anciaux <guillaume.anciaux@epfl.ch>",
]
__copyright__ = "Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__license__ = "LGPLv3"

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    prank = comm.Get_rank()
except ImportError:
    prank = 0

import akantu as aka
import numpy as np


# -----------------------------------------------------------------------------
def solve(material_file, mesh_file, body_stress):
    aka.parseInput(material_file)
    spatial_dimension = 2

    # -------------------------------------------------------------------------
    # Initialization
    # -------------------------------------------------------------------------
    mesh = aka.Mesh(spatial_dimension)
    if prank == 0:
        mesh.read(mesh_file)

    mesh.distribute()

    model = aka.SolidMechanicsModel(mesh)
    model.initFull(_analysis_method=aka._static)

    model.setBaseName("square")
    model.addDumpFieldVector("displacement")
    model.addDumpFieldVector("external_force")
    model.addDumpField("strain")
    model.addDumpField("stress")
    model.addDumpField("blocked_dofs")

    # -------------------------------------------------------------------------
    # Boundary conditions
    # -------------------------------------------------------------------------
    model.applyBC(aka.FixedValue(0.0, aka._x), "Fixed_y")
    model.applyBC(aka.FixedValue(0.0, aka._y), "Fixed_y")

    gravity = np.eye(2)*body_stress
    model.applyBC(aka.BodyStressFunctor(gravity), "steel")
    
    solver = model.getNonLinearSolver()
    solver.set("max_iterations", 2)
    solver.set("threshold", 1e-10)
    solver.set("convergence_type", aka.SolveConvergenceCriteria.solution)

    model.solveStep()

    model.dump()


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def main():
    mesh_file = 'square.msh'
    material_file = 'material.dat'
    body_stress = 1.

    solve(material_file, mesh_file, body_stress)


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
