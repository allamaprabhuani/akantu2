#!/usr/bin/env python3

""" custom_constitutive_law.py: Custom constitutive law example"""

__author__ = "Mohit Pundir"
__credits__ = [
    "Mohit Pundir <mpundir@ethz.ch>",
]
__copyright__ = "Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__license__ = "LGPLv3"

import numpy as np
import akantu as aka


class LocalConstitutiveLaw(aka.ConstitutiveLaw):

    def __init__(self, model, _id):
        super().__init__(model, _id)
        super().registerParamReal('permitivity', aka._pat_readable |
                                  aka._pat_parsable, 'permitivity')

    def initConstitutiveLaw(self):
        super().initConstitutiveLaw()

    def getCelerity(self):
        return 8.* self.getReal('permitivity')

    def computeFlux(self, el_type, ghost_type):
        concentration = self.getModel().getDof()
        concentration_flux = self.getFluxDof(el_type, ghost_type)
        concentration_gradient = self.getGradientDof(el_type, ghost_type)

        fem = self.getModel().getFEEngine()
        fem.gradientOnIntegrationPoints(concentration, concentration_gradient,
                                        1, el_type, ghost_type)
        concentration_flux[:, :, :] = np.einsum('a,ij->aij', self.getReal('D'),
                                                concentration_gradient)



    def computeTangentModulii(self, el_type, tangent_matrix, ghost_type):
        n_quads = tangent_matrix.shape[0]
        tangent = tangent_matrix.reshape(n_quads, 1, 1)

        tangent[:, 0, 0] = self.getReal('permitivity')

        print(tangent)

    def getEffectiveCapacity(self):
        return 1.


# register constitutive law to the ConstitutiveLawFactory
def allocator(unused, model, _id):
    return LocalConstitutiveLaw(model, _id)



def solve(material_file, mesh_file):
    aka.parseInput(material_file)
    spatial_dimension = 1

    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = aka.PoissonModel(mesh)
    model.initFull()

    time_step = model.getStableTimeStep() * 0.1
    model.setTimeStep(time_step)

    model.setBaseName('custom')
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
    law_factory = aka.ConstitutiveLawFactory.getInstance()
    law_factory.registerAllocator("local_constitutive_law", allocator)


    mesh_file = 'bar.msh'
    material_file = 'custom_law.dat'

    solve(material_file, mesh_file)


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
