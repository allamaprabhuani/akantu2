/**
 * @file   material_viscoplasticity.cc
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 *
 * @date   Tue Jul 09 18:15:37 20130
 *
 * @brief  Implementation of the inline functions of the material viscoplasticity
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
#include "material_viscoplasticity.hh"


/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void MaterialViscoPlasticity<dim>::computeStressOnQuad(Matrix<Real> & grad_u, Matrix<Real> & grad_delta_u, Matrix<Real> & sigma, Matrix<Real> & inelas_strain, Real & iso_hardening) {
    //Infinitesimal plasticity

  Matrix<Real> d_inelas_strain(dim,dim);
  //Real r=iso_hardening;
  
  // Compute stress magnitude

  Real s = sigma.doubleDot(sigma);
  Real sigma_mag=sqrt(s);

  //Compute plastic strain increment

  for (UInt i = 0; i < dim; ++i) 
    for (UInt j = 0; j < dim; ++j) 
      d_inelas_strain(i,j) = (ts * edot0 * pow(sigma_mag,(rate-1)) / pow(sigmay+iso_hardening,rate)) * sigma(i,j);  

 // Compute plastic strain increment magnitude

  s = d_inelas_strain.doubleDot(d_inelas_strain);
  Real dep_mag=sqrt(s);

  //Update stress and plastic strain

  for (UInt i = 0; i < dim; ++i) 
    for (UInt j = 0; j < dim; ++j) {    
      sigma(i,j) = sigma(i,j) + 2. * mu * (0.5 * (grad_delta_u(i,j)+ grad_delta_u(j,i))-d_inelas_strain(i,j)) + (i==j) * lambda * (grad_delta_u.trace() - d_inelas_strain.trace());
      inelas_strain(i,j) = inelas_strain(i,j) +  d_inelas_strain(i,j);
    }

  //Update resistance stress
  
  iso_hardening = iso_hardening + h * dep_mag;

} 


/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void MaterialViscoPlasticity<dim>::computeTangentModuliOnQuad(Matrix<Real> & tangent, Matrix<Real> & grad_delta_u, Matrix<Real> & sigma_tensor, Matrix<Real>  & previous_sigma_tensor, Real & iso_hardening) {

  UInt cols = tangent.cols();
  UInt rows = tangent.rows();

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

	      tangent(m,n) = (i==k) * (j==l) * 2. * mu + (i==j) * (k==l) * lambda;
	      if ((m==n) && (m>=dim))
		tangent(m,n) = tangent(m,n) - mu;

	}
 }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline Real MaterialViscoPlasticity<spatial_dimension>::getStableTimeStep(Real h,
        __attribute__((unused)) const Element & element) {
    return (h / getPushWaveSpeed());
}
