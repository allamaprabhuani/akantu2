#!/usr/bin/env python3

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
