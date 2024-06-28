#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import sys
sys.path.insert(0,"/home/ble/developpement/akantu/build/buildCZMDamageDebug/python")
print(sys.path)

import akantu as aka
import numpy as np
import matplotlib.pyplot as plt

from czm_damage import *
					
def set_dumpers(model):
    model.setBaseName("cohesive")
    model.addDumpFieldVector("displacement")
    model.addDumpFieldVector("external_force")
    model.addDumpField("strain")
    model.addDumpField("stress")
    model.addDumpField("blocked_dofs")

    model.setBaseNameToDumper("cohesive elements", "cohesive")
    model.addDumpFieldVectorToDumper("cohesive elements", "displacement")
    model.addDumpFieldToDumper("cohesive elements", "damage")
    model.addDumpFieldVectorToDumper("cohesive elements", "tractions")
    model.addDumpFieldVectorToDumper("cohesive elements", "opening")


class FixedDisplacement (aka.DirichletFunctor):
    '''
        Fix the displacement at its current value
    '''

    def __init__(self, axis, vel):
        super().__init__(axis)
        self.axis = axis
        self.time = 0
        self.vel = vel

    def set_time(self, t):
        self.time = t

    def get_imposed_disp(self):
        return self.vel*self.time
        
    def __call__(self, node, flags, disp, coord):
        # sets the blocked dofs vector to true in the desired axis
        flags[int(self.axis)] = True
        disp[int(self.axis)] = self.get_imposed_disp()
        
def solve(material_file, mesh_file):
    aka.parseInput(material_file)
    spatial_dimension = 2
    L = 0.4
    
    # -------------------------------------------------------------------------
    # Initialization
    # -------------------------------------------------------------------------
    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = aka.SolidMechanicsModelCohesive(mesh)
    model.getElementInserter().setLimit(aka._x, 0.5*L-0.1, 0.5*L+0.1);

    model.initFull(_analysis_method=aka._static, _is_extrinsic=True)
    
    # Initilize a new solver (explicit Newmark with lumped mass)
    model.initNewSolver(aka._explicit_lumped_mass)
    # Dynamic insertion of cohesive elements
    model.updateAutomaticInsertion()
                       
    set_dumpers(model)


    E = model.getMaterial("steel").getReal("E")
    Gc = model.getMaterial("cohesive").getReal("G_c")
    sigc = model.getMaterial("cohesive").getReal("sigma_c")
    wc = 2.*Gc/sigc
    w0 = (sigc/E)*L
    k = model.getMaterial("cohesive").getReal("k")
    g = degradation_function_linear_czm()
    h = softening_function_linear_czm(k,sigc,Gc)
    updater = damage_updater(Gc,sigc,k,g,h)

    imposed_disp = [0.]
    reaction_force = [0.]
    E_pot = []
    E_kin = []
    E_dis = []
    E_rev = []
    E_con = []

    U = 0.
    F = 0.
        
    # -------------------------------------------------------------------------
    # Initial and boundary conditions
    # -------------------------------------------------------------------------
    eps0dot = 5e-1
    tc = w0/(eps0dot*L)
    functor_r = FixedDisplacement(aka._x, L*eps0dot)
    functor_r.set_time(tc)
    model.applyBC(functor_r, 'right')

    model.applyBC(aka.FixedValue(0.0, aka._x), "left")
    model.applyBC(aka.FixedValue(0.0, aka._y), "point")
    
    nodes = model.getMesh().getNodes()
    disp_field = np.zeros(nodes.shape)
    disp_field[:, 0] = nodes[:, 0]*(w0/L)
    model.getDisplacement()[:] = disp_field
    vel_field = np.zeros(nodes.shape)
    vel_field[:, 0] = nodes[:, 0]*eps0dot
    model.getVelocity()[:] = vel_field
    
    model.getExternalForce()[:] = 0
    model.assembleInternalForces()                           
                                                  
    d = model.getMaterial("cohesive").getInternalReal("czm_damage")(aka._segment_2)
    lda = model.getMaterial("cohesive").getInternalReal("lambda")(aka._segment_2)
    Fint = model.getInternalForce() 

    nodes_right = mesh.getElementGroup("right").getNodeGroup().getNodes()    
    
    U = functor_r.get_imposed_disp()
    imposed_disp.append(U)    
    F = -np.sum(Fint[nodes_right,0])
    reaction_force.append(F)
    
    Ep = model.getEnergy("potential")
    Ek = model.getEnergy("kinetic")
    Ed = model.getEnergy("dissipated")
    Er = model.getEnergy("reversible")
    
    E_pot.append(Ep)
    E_kin.append(Ek)
    E_dis.append(Ed)
    E_rev.append(Er)
    
    model.dump()
    #model.dump("cohesive elements")
            
    maxsteps = 2
    dt = model.getStableTimeStep()*0.1
    # choose the timestep
    model.setTimeStep(dt)

    time = tc
    for i in range(0, maxsteps):
        print("---- step {0}/{1}".format(i, maxsteps))
        time+=dt
        functor_r.set_time(time)
        # fix displacements of the right boundary
        model.applyBC(functor_r, 'right')
        d_previous = model.getMaterial("cohesive").getInternalReal("czm_damage")(aka._segment_2).copy()        
        #print("lda = ",lda)
        model.solveStep('explicit_lumped')
        model.computeLagrangeMultiplier()
        d[:] = updater.update_d(d_previous,lda)
        #print("d = ",d)

        model.checkCohesiveInsertion()
        
        U = functor_r.get_imposed_disp()
        imposed_disp.append(U)    
        F = -np.sum(Fint[nodes_right,0])
        reaction_force.append(F)

        Ep = model.getEnergy("potential")
        Ek = model.getEnergy("kinetic")
        Ed = model.getEnergy("dissipated")
        Er = model.getEnergy("reversible")
    
        E_pot.append(Ep)
        E_kin.append(Ek)
        E_dis.append(Ed)
        E_rev.append(Er)
        
        model.dump()
        #model.dump("cohesive elements")

    plt.plot(imposed_disp,reaction_force,'-o')
    plt.plot([0., w0, wc], [0., sigc*L, 0.], 'k-o')
    plt.show()


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def main():
    mesh_file = "triangle.msh"
    material_file = "material.dat"
    solve(material_file, mesh_file)


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
