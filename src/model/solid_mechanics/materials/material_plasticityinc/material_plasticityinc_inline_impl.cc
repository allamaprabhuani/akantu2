/**
 * @file   material_plasticity_inc.cc
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 *
 * @date   Tue Jul 09 18:15:37 20130
 *
 * @brief  Implementation of the inline functions of the material plasticity
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <cmath>
#include "material_plasticityinc.hh"


/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void MaterialPlasticityinc<dim>::computeStressOnQuad(Matrix<Real> & grad_u, Matrix<Real> & grad_delta_u, Matrix<Real> & sigma, Matrix<Real> & inelas_strain, Real & iso_hardening) {
    //Infinitesimal plasticity
  Matrix<Real> sigma_tr(dim, dim);
  Matrix<Real> sigma_tr_dev(dim, dim);
  Matrix<Real> d_inelas_strain(dim,dim);
  //Real r=iso_hardening;
  Real dp=0.0;
  Real d_dp=0.0;
  UInt n=0;

  //Compute trial stress, sigma_tr

  for (UInt i = 0; i < dim; ++i)
        for (UInt j = 0; j < dim; ++j)
	  sigma_tr(i,j) = sigma(i,j) + 2. * mu * 0.50*(grad_delta_u(i,j)+grad_delta_u(j,i)) + (i == j) * lambda * grad_delta_u.trace();

  // Compute deviatoric trial stress,  sigma_tr_dev

for (UInt i = 0; i < dim; ++i)
        for (UInt j = 0; j < dim; ++j)
	  sigma_tr_dev(i,j) = sigma_tr(i,j) - (i == j) * sigma_tr.trace() / 3.0;

// Compute effective deviatoric trial stress

  Real s = sigma_tr_dev.doubleDot(sigma_tr_dev);
  Real sigma_tr_dev_eff=sqrt(3./2. * s);


   //Loop for correcting stress based on yield function

  while ((sigma_tr_dev_eff-iso_hardening-sigmay) > 0){
    //r = r +  h * dp;

    d_dp = (sigma_tr_dev_eff - 3. * mu *dp -  iso_hardening - sigmay) / (3. * mu + h);
    iso_hardening = iso_hardening + h * d_dp;
    dp = dp + d_dp;
    n = n+1;

          
    if ((d_dp < 1e-5)|| (n>50))
      break;
  }

 //Update internal variable

    for (UInt i = 0; i < dim; ++i) {
      for (UInt j = 0; j < dim; ++j) {
     
	if (sigma_tr_dev_eff)

	  d_inelas_strain(i,j) = 3. / 2. * dp * sigma_tr_dev(i,j) / sigma_tr_dev_eff;
	  inelas_strain(i,j) = inelas_strain(i,j) +  d_inelas_strain(i,j);
	  sigma(i,j) = sigma(i,j) + 2. * mu * (0.5 * (grad_delta_u(i,j)+ grad_delta_u(j,i))-d_inelas_strain(i,j)) + (i==j) * lambda * (grad_delta_u.trace() - d_inelas_strain.trace());

      }
    }
} 


/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void MaterialPlasticityinc<dim>::computeTangentModuliOnQuad(Matrix<Real> & tangent, Matrix<Real> & grad_delta_u, Matrix<Real> & sigma_tensor, Matrix<Real>  & previous_sigma_tensor, Real & iso_hardening) {

  UInt cols = tangent.cols();
  UInt rows = tangent.rows();
  Matrix<Real> sigma_dev(dim, dim);
  Matrix<Real> sigma_tr(dim, dim);
  Matrix<Real> sigma_tr_dev(dim, dim);
  // Real r=iso_hardening;
 

  //Compute trial stress, sigma_tr

  for (UInt i = 0; i < dim; ++i)
        for (UInt j = 0; j < dim; ++j)
	  	  sigma_tr(i,j) = previous_sigma_tensor(i,j) + 2. * mu * 0.50*(grad_delta_u(i,j)+grad_delta_u(j,i)) + (i == j) * lambda * grad_delta_u.trace();

  // Compute deviatoric trial stress,  sigma_tr_dev

for (UInt i = 0; i < dim; ++i)
        for (UInt j = 0; j < dim; ++j)
	  sigma_tr_dev(i,j) = sigma_tr(i,j) - (i == j) * sigma_tr.trace() / 3.0;

// Compute effective deviatoric trial stress

  Real s = sigma_tr_dev.doubleDot(sigma_tr_dev);
  Real sigma_tr_dev_eff=sqrt(3./2. * s);

  // Compute deviatoric stress,  sigma_dev

for (UInt i = 0; i < dim; ++i)
        for (UInt j = 0; j < dim; ++j)
       	  sigma_dev(i,j) = sigma_tensor(i,j) - (i == j) * sigma_tensor.trace() / 3.0;

// Compute effective deviatoric stress

  s =  sigma_dev.doubleDot(sigma_dev);
  Real sigma_dev_eff=sqrt(3./2. * s);

  Real xr = 0.0;
  if(sigma_tr_dev_eff) 
    xr = sigma_dev_eff / sigma_tr_dev_eff ;

  Real q = 1.5 * (1. / (1. +  3. * mu  / h) - xr);

for (UInt m = 0; m < rows; m++) {
  UInt i, j;
        if (m < dim) {
            i = m;
            j = m;
        } else {
            if (dim == 3) {
                if (m == 3) {
                    i = 1;
                    j = 2;
                } else if (m == 4) {
                    i = 0;
                    j = 2;
                } else if (m == 5) {
                    i = 0;
                    j = 1;
                }
            } else if (dim == 2) {
                if (m == 2) {
                    i = 0;
                    j = 1;
                }
            }
        }

        for (UInt n = 0; n < cols; n++) {
            UInt k, l;
            if (n < dim) {
                k = n;
                l = n;
            } else {
                if (dim == 3) {
                    if (n == 3) {
                        k = 1;
                        l = 2;
                    } else if (n == 4) {
                        k = 0;
                        l = 2;
                    } else if (n == 5) {
                        k = 0;
                        l = 1;
                    }
                } else if (dim == 2) {
                    if (n == 2) {
                        k = 0;
                        l = 1;
                    }
                }
            }

	    	   if (((sigma_tr_dev_eff-iso_hardening-sigmay) > 0) && (xr > 0)){

	       tangent(m,n) = 2. * mu * q * (sigma_tr_dev (i,j) / sigma_tr_dev_eff) * (sigma_tr_dev (k,l) / sigma_tr_dev_eff) + (i==k) * (j==l) * 2. * mu * xr + (i==j) * (k==l)* (kpa - 2./3. * mu * xr); 
	       if ((m==n) && (m>=dim))
		 tangent(m,n) = tangent(m,n) - mu * xr;

	    }

	    else {

	      tangent(m,n) = (i==k) * (j==l) * 2. * mu + (i==j) * (k==l) * lambda;
	      if ((m==n) && (m>=dim))
		tangent(m,n) = tangent(m,n) - mu;

	        }
 

	    //correct tangent stiffness for shear component



        }
 }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline Real MaterialPlasticityinc<spatial_dimension>::getStableTimeStep(Real h,
        __attribute__((unused)) const Element & element) {
    return (h / getPushWaveSpeed());
}
