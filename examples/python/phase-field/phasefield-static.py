#!/usr/bin/env python
# coding: utf-8
""" phasefield-static.py: Static phase field example"""

__author__ = "Mohit Pundir"
__credits__ = [
    "Mohit Pundir <mohit.pundir@epfl.ch>",
]
__copyright__ = "Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__license__ = "LGPLv3"


import numpy as np
import akantu as aka


aka.parseInput("material_static.dat")

dim = 2
mesh = aka.Mesh(dim)
mesh.read("plate_static.msh")

model = aka.CouplerSolidPhaseField(mesh)

solid = model.getSolidMechanicsModel()
phase = model.getPhaseFieldModel()

solid.initFull(_analysis_method=aka._static)
solver = solid.getNonLinearSolver('static')
solver.set('max_iterations', 100)
solver.set('threshold', 1e-9)
solver.set("convergence_type", aka.SolveConvergenceCriteria.solution)


solid.getNewSolver("linear_static", aka.TimeStepSolverType.static,
                   aka.NonLinearSolverType.linear)
solid.setIntegrationScheme("linear_static", "displacement",
                           aka.IntegrationSchemeType.pseudo_time)


phase.initFull(_analysis_method=aka._static)
phase.getNewSolver("nonlinear_static", aka.TimeStepSolverType.static,
                   aka.NonLinearSolverType.newton_raphson)
phase.setIntegrationScheme("nonlinear_static", "damage",
                           aka.IntegrationSchemeType.pseudo_time)

solver = phase.getNonLinearSolver('nonlinear_static')
solver.set('max_iterations', 100)
solver.set('threshold', 1e-4)
solver.set("convergence_type", aka.SolveConvergenceCriteria.solution)


solid.applyBC(aka.FixedValue(0, aka._y), "bottom")
solid.applyBC(aka.FixedValue(0, aka._x), "left")

# Initialization for bulk vizualisation
solid.setBaseName('phasefield-static')
solid.addDumpFieldVector('displacement')
solid.addDumpFieldVector('external_force')
solid.addDumpField('strain')
solid.addDumpField('stress')
solid.addDumpField('damage')
solid.addDumpField('blocked_dofs')

nb_dofs = solid.getMesh().getNbNodes() * dim

increment = solid.getIncrement()
displacement = solid.getDisplacement()
displacement = displacement.reshape(nb_dofs)

blocked_dofs = solid.getBlockedDOFs()
blocked_dofs = blocked_dofs.reshape(nb_dofs)

damage = phase.getDamage()

tolerance = 1e-6

steps = 1500
increment = 1e-5

for n in range(steps):
    print("Computing iteration " + str(n + 1) + "/" + str(steps))

    solid.applyBC(aka.IncrementValue(increment, aka._y), 'top')

    mask = blocked_dofs == False  # NOQA: E712

    iiter = 0
    error_disp = 1
    error_dam = 1

    displacement_prev = displacement[mask].copy()

    damage_prev = damage.copy()
    damage_prev = damage_prev

    # solve using staggered scheme
    while (error_disp > tolerance or error_dam > tolerance):
        model.solve("linear_static", "")

        displacement_new = displacement[mask]
        damage_new = damage

        delta_disp = displacement_new - displacement_prev
        delta_dam = damage_new - damage_prev

        error_disp = np.linalg.norm(delta_disp)
        error_dam = np.linalg.norm(delta_dam)

        iiter += 1

        displacement_prev = displacement_new.copy()
        damage_prev = damage_new.copy()

        print(error_dam, error_disp)
        if iiter > 500:
            raise Exception('Convergence not reached')

    if n % 50 == 0:
        solid.dump()

solid.dump()
