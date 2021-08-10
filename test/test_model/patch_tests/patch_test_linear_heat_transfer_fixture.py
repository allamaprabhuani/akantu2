#!/usr/bin/env python3

"""patch_test_linear_heat_transfer_fixture.py: heat transfer patch test in
python"""

__author__ = "Guillaume Anciaux"
__credits__ = [
    "Guillaume Anciaux <guillaume.anciaux@epfl.ch>",
]
__copyright__ = "Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__license__ = "LGPLv3"

import patch_test_linear_fixture
import akantu


class TestPatchTestHTMLinear(patch_test_linear_fixture.TestPatchTestLinear):

    model_type = akantu.HeatTransferModel

    def applyBC(self):
        super().applyBC()
        temperature = self.model.getTemperature()
        self.applyBConDOFs(temperature)

    def checkAll(self):
        temperature = self.model.getTemperature()
        C = self.model.getMatrix("conductivity")
        self.checkDOFs(temperature)
        self.checkGradient(self.model.getTemperatureGradient(self.elem_type),
                           temperature)

        self.prescribed_gradient(temperature)
        self.checkResults(lambda grad_T: C.dot(grad_T.T),
                          self.model.getKgradT(self.elem_type),
                          temperature)

    def initModel(self, method, material_file):
        super().initModel(method, material_file)

        if method != akantu._static:
            self.model.setTimeStep(0.5 * self.model.getStableTimeStep())


def run_test_generic(self_, method):
    self_.initModel(method, "heat_transfer_input.dat")

    coordinates = self_.mesh.getNodes()
    temperature = self_.model.getTemperature()
    #  set the position of all nodes to the static solution
    self_.setLinearDOF(temperature, coordinates)

    for s in range(0, 100):
        self_.model.solveStep()

    self_.checkAll()
