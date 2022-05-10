#!/usr/bin/env python3
""" diffusion_law.py: Poisson model with diffusion example"""

__author__ = "Mohit Pundir"
__credits__ = [
    "Mohit Pundir <mpundir@ethz.ch>",
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


def solve(material_file, mesh_file):
    aka.parseInput(material_file)
    spatial_dimension = 1

    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = aka.PoissonModel(mesh)
    model.initFull()

    time_step = model.getStableTimeStep() * 0.1
    model.setTimeStep(time_step)

    model.setBaseName('diffusion')
    model.addDumpField('dof')
    model.addDumpField('internal_dof_rate')
    model.addDumpField('external_dof_rate')
    model.dump()

    external_flux_rate = model.getExternalDofRate()
    external_flux_rate[0] = 1e-8

    nb_steps = 1000
    for s in range(nb_steps):
        model.solveStep()

        if s % 100 == 0:
            model.dump()

        print('Step {0}/{1}'.format(s, nb_steps))


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def main():
    mesh_file = 'bar.msh'
    material_file = 'diffusion.dat'

    solve(material_file, mesh_file)


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
