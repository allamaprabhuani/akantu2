#!/usr/bin/env python3

""" test_patch_linear_elastic_static.py: static patch test in python"""

__author__ = "Guillaume Anciaux"
__credits__ = [
    "Guillaume Anciaux <guillaume.anciaux@epfl.ch>",
]
__copyright__ = "Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__license__ = "LGPLv3"

import sys
from patch_test_linear_solid_mechanics_fixture import TestPatchTestSMMLinear
import akantu


def foo(self):
    filename = "material_check_stress_plane_stress.dat"
    if self.plane_strain:
        filename = "material_check_stress_plane_strain.dat"

    self.initModel(akantu._static, filename)

    solver = self.model.getNonLinearSolver()
    solver.set("max_iterations", 2)
    solver.set("threshold", 2e-4)
    solver.set("convergence_type", akantu.SolveConvergenceCriteria.residual)

    self.model.solveStep()

    self.checkAll()


def test():
    TestPatchTestSMMLinear.TYPED_TEST(foo, "Static")


if 'pytest' not in sys.modules:
    test()
