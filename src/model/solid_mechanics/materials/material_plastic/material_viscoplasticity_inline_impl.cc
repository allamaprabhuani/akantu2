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
inline void MaterialViscoPlasticity<dim>::computeStressOnQuad(Matrix<Real> & grad_u,
                                                              Matrix<Real> & grad_delta_u,
                                                              Matrix<Real> & sigma,
                                                              Matrix<Real> & inelas_strain,
                                                              Real & iso_hardening) {
  // Infinitesimal plasticity
  Matrix<Real> d_inelas_strain(dim,dim);

  // Compute stress magnitude
  Real s = sigma.doubleDot(sigma);
  Real sigma_mag=sqrt(s);

  // Compute plastic strain increment
  for (UInt i = 0; i < dim; ++i) {
    for (UInt j = 0; j < dim; ++j) {
      d_inelas_strain(i,j) = (ts * edot0 * pow(sigma_mag,(rate-1)) / pow(sigmay+iso_hardening,rate)) * sigma(i,j);
    }
  }

  // Compute plastic strain increment magnitude
  s = d_inelas_strain.doubleDot(d_inelas_strain);
  Real dep_mag = std::sqrt(s);

  // Update stress and plastic strain
  Matrix<Real> grad_u_elsatic(dim, dim);
  grad_u_elsatic = grad_delta_u;
  grad_u_elsatic -= d_inelas_strain;

  Matrix<Real> sigma_elastic(dim, dim);
  MaterialElastic<dim>::computeStressOnQuad(grad_u_elsatic, sigma_elastic);
  sigma += sigma_elastic;

  inelas_strain += d_inelas_strain;

  //Update resistance stress
  iso_hardening = iso_hardening + h * dep_mag;
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void MaterialViscoPlasticity<dim>::computeTangentModuliOnQuad(Matrix<Real> & tangent,
                                                                     Matrix<Real> & grad_delta_u,
                                                                     Matrix<Real> & sigma_tensor,
                                                                     Matrix<Real> & previous_sigma_tensor,
                                                                     Real & iso_hardening) {
  UInt cols = tangent.cols();
  UInt rows = tangent.rows();

  for (UInt m = 0; m < rows; ++m) {
    UInt i = VoigtHelper<dim>::vec[m][0];
    UInt j = VoigtHelper<dim>::vec[m][1];

    for (UInt n = 0; n < cols; ++n) {
      UInt k = VoigtHelper<dim>::vec[n][0];
      UInt l = VoigtHelper<dim>::vec[n][1];
      tangent(m,n) = (i==k) * (j==l) * 2. * this->mu + (i==j) * (k==l) * this->lambda;
      tangent(m,n) -= (m==n) * (m>=dim) * this->mu;
    }
  }
}
