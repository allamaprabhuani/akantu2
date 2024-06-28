#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import sys
sys.path.insert(0,"/home/ble/buildAkantu/buildDebug/python")
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


def solve(material_file, mesh_file, increment):
    aka.parseInput(material_file)
    spatial_dimension = 2

    # -------------------------------------------------------------------------
    # Initialization
    # -------------------------------------------------------------------------
    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = aka.SolidMechanicsModelCohesive(mesh)
    model.getElementInserter().setLimit(aka._x, -0.1, 0.1);

    model.initFull(_analysis_method=aka._static, _is_extrinsic=False)
                           
    set_dumpers(model)

    L = 0.4
    E = model.getMaterial("steel").getReal("E")
    Gc = model.getMaterial("cohesive").getReal("G_c")
    tc = model.getMaterial("cohesive").getReal("sigma_c")
    wc = 2.*Gc/tc
    k = model.getMaterial("cohesive").getReal("k")
    g = degradation_function_linear_czm()
    h = softening_function_linear_czm(k,tc,Gc)
    updater = damage_updater(Gc,tc,k,g,h)

    imposed_disp = [0.]
    reaction_force = [0.]
    U = 0.
    F = 0.
        
    # -------------------------------------------------------------------------
    # Boundary conditions
    # -------------------------------------------------------------------------
    model.applyBC(aka.FixedValue(0.0, aka._x), "left")
    model.applyBC(aka.FixedValue(0.0, aka._y), "point")
    model.applyBC(aka.IncrementValue(increment, aka._x), "right")  
    U+=increment
    imposed_disp.append(U)
    
    model.getExternalForce()[:] = 0

    solver = model.getNonLinearSolver("static")
    solver.set("max_iterations", 2)
    solver.set("threshold", 1e-10)
    solver.set("convergence_type", aka.SolveConvergenceCriteria.residual)

    model.getNewSolver("linear_static", aka.TimeStepSolverType.static,
                       aka.NonLinearSolverType.linear)
    model.setIntegrationScheme("linear_static", "displacement",
                       aka.IntegrationSchemeType.pseudo_time)
    model.setIntegrationScheme("linear_static", "lambda",
                       aka.IntegrationSchemeType.pseudo_time)
                                                  
    d = model.getMaterial("cohesive").getInternalReal("czm_damage")(aka._cohesive_2d_4)
    lda = model.getMaterial("cohesive").getInternalReal("lambda")(aka._cohesive_2d_4)
    Fint = model.getInternalForce() 
    
    nodes_right = mesh.getElementGroup("right").getNodeGroup().getNodes()
    model.solveStep("linear_static")
    F = -np.sum(Fint[nodes_right,0])
    reaction_force.append(F)
    
    model.dump()
    #model.dump("cohesive elements")
            
    maxsteps = 8
    maxiter = 10
    tol = 1e-6
    
    for i in range(0, maxsteps):
        print("---- step {0}/{1}".format(i, maxsteps))
        model.applyBC(aka.IncrementValue(increment, aka._x), "right")
        U+=increment
        imposed_disp.append(U)    
        it = 0
        d_previous = model.getMaterial("cohesive").getInternalReal("czm_damage")(aka._cohesive_2d_4).copy()
        W0 = energy(model)
        while it < maxiter:
            #aka.debug.setDebugLevel(aka.dblInfo);
            model.solveStep("linear_static")
            d[:] = updater.update_d(d_previous,lda)
            W = energy(model)
            err = abs(W-W0)/W
            print("iter ",it,", err = ",err)          
            if err < tol:
                break
            else:
                W0 = W
                it+=1
        F = -np.sum(Fint[nodes_right,0])
        reaction_force.append(F)
    
        model.dump()
        #model.dump("cohesive elements")

    plt.plot(imposed_disp,reaction_force,'-o')
    plt.plot([0., (tc/E)*L, wc], [0., tc*L, 0.], 'k-o')
    plt.show()
    W = model.getEnergy("potential");
    D = model.getEnergy("dissipated");
    print("Elastic energy = ",W)
    print("Dissipated energy = ",D)

# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def main():
    mesh_file = "triangle.msh"
    material_file = "material.dat"
    increment = 0.01e-3
    solve(material_file, mesh_file, increment)


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
