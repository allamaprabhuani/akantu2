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

import scipy.optimize

import matplotlib.pyplot as plt

class degradation_function:
	# g(d)
	def val(self,d):
		return
	# dg_dd(d)
	def deriv(self,d):
		return
	# 1+g(d)
	def augmented_stiffness(self,d):
		return
	def augmented_compliance(self,d):
	# 1/(1+g(d))
		return
	def augmented_compliance_deriv(self,d):
	# d()1/(1+g(d))_dd
		return
				
class degradation_function_linear_czm(degradation_function):
	def val(self,d):
		return 1./d-1.
	def deriv(self,d):
		return -1./d**2
	def augmented_stiffness(self,d):
		return 1./d
	def augmented_compliance(self,d):
		return d
	def augmented_compliance_deriv(self,d):
		return 1.
		
class softening_function:
	def val(self,d):
		return
	def deriv(self,d):
		return
		
class softening_function_linear_czm:
	def __init__(self,k,tc,Gc):
		wc = 2.*Gc/tc
		self.a = k*wc/tc
	def val(self,d):
		return d/(d+self.a*(1.-d))
	def deriv(self,d):
		return self.a/(d+self.a*(1.-d))**2
		
class damage_updater:
    def __init__(self,Gc, sigc, k, g, h):
        self.Gc = Gc
        self.sigc = sigc
        self.k = k
        self.g = g
        self.h = h
    def update_d(self,d_previous,lda):
        dz = d_previous.flatten()
        shape_d = d_previous.shape
        lda_lda_on_k = (lda[:,0]**2 + lda[:,1]**2)/self.k        
        
        def F(d):		
            return -0.5*np.dot(lda_lda_on_k,self.g.augmented_compliance(d))+self.Gc*np.sum(self.h.val(d))
        def dF(d):
            return -0.5*self.g.augmented_compliance_deriv(d)*lda_lda_on_k+self.Gc*self.h.deriv(d)

        #dd = np.linspace(0.,1.,100)
        #D = dz.copy()
        #FF = []
        #for i in range(len(dd)):
        #    D[0] = dd[i]
        #    FF.append(F(D))
        #plt.plot(dd,FF)
        #plt.show()
        #print("dz = ",dz)
        #print("lda_lda_on_k = ",lda_lda_on_k)
        #print("F1  = ",-0.5*np.dot(lda_lda_on_k,self.g.augmented_compliance(dz)))
        #print("F2 = ",self.Gc*np.sum(self.h.val(dz)))
        #print("F(dz) = ",F(dz))
        #print("self.k = ",self.k)
        
        boundd = scipy.optimize.Bounds(dz, np.ones(len(dz))) 
        options = {'verbose':1}
        resd = scipy.optimize.minimize(F, dz, bounds=boundd
		                              , jac =dF
                                      #, method = 'trust-constr'
                                      #, options = options
                                      #, tol = 1e-10
                                      )
        if not resd.success:
            print("resd.message = ",resd.message)
        #print("niterd = ",resd.nit)
        return np.resize(resd.x,shape_d)
	
def energy(model):
    return model.getEnergy("potential") + model.getEnergy("dissipated")
					
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
