/**
 * @file   material_plasticity_inc.cc
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
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
inline void MaterialPlasticityinc<dim>::computeStressOnQuad(Matrix<Real> & grad_u,
							    Matrix<Real> & grad_delta_u,
							    Matrix<Real> & sigma,
							    Matrix<Real> & inelas_strain,
							    Real & iso_hardening,
							    Real sigma_th_cur,
							    Real sigma_th_prev) {
  //Infinitesimal plasticity
  Matrix<Real> sigma_tr(dim, dim);
  Matrix<Real> sigma_tr_dev(dim, dim);
  Matrix<Real> d_inelas_strain(dim,dim);
  //Real r=iso_hardening;
  Real dp=0.0;
  Real d_dp=0.0;
  UInt n=0;

  Real delta_sigma_th = sigma_th_cur - sigma_th_prev;

  //Compute trial stress, sigma_tr
  MaterialElastic<dim>::computeStressOnQuad(sigma, grad_delta_u, delta_sigma_th);
  sigma_tr += sigma;

  // Compute deviatoric trial stress,  sigma_tr_dev
  sigma_tr_dev  = sigma_tr;
  sigma_tr_dev -= Matrix<Real>::eye(dim, sigma_tr.trace() / 3.0);

  // Compute effective deviatoric trial stress
  Real s = sigma_tr_dev.doubleDot(sigma_tr_dev);
  Real sigma_tr_dev_eff = std::sqrt(3./2. * s);

  //Loop for correcting stress based on yield function
  while ((sigma_tr_dev_eff-iso_hardening-sigmay) > 0) {
    //r = r +  h * dp;
    d_dp = (sigma_tr_dev_eff - 3. * this->mu *dp -  iso_hardening - sigmay) / (3. * this->mu + h);

    iso_hardening = iso_hardening + h * d_dp;
    dp = dp + d_dp;

    ++n;
    if ((d_dp < 1e-5) || (n>50))
      break;
  }

  //Update internal variable
  if (std::abs(sigma_tr_dev_eff) >
      sigma_tr_dev.norm<L_inf>() * Math::getTolerance()) {
    d_inelas_strain = sigma_tr_dev;
    d_inelas_strain *= 3./2. * dp / sigma_tr_dev_eff;
  }

  Matrix<Real> grad_u_elsatic(dim, dim);
  grad_u_elsatic  = grad_delta_u;
  grad_u_elsatic -= d_inelas_strain;

  Matrix<Real> sigma_elastic(dim, dim);
  MaterialElastic<dim>::computeStressOnQuad(sigma_elastic, grad_u_elsatic, delta_sigma_th);
  sigma += sigma_elastic;

  inelas_strain += d_inelas_strain;
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void MaterialPlasticityinc<dim>::computeTangentModuliOnQuad(Matrix<Real> & tangent,
								   Matrix<Real> & grad_delta_u,
								   Matrix<Real> & sigma_tensor,
								   Matrix<Real> & previous_sigma_tensor,
								   Real & iso_hardening) {
  UInt cols = tangent.cols();
  UInt rows = tangent.rows();
  Matrix<Real> sigma_dev(dim, dim);
  Matrix<Real> sigma_tr(dim, dim);
  Matrix<Real> sigma_tr_dev(dim, dim);
  // Real r=iso_hardening;

  //Compute trial stress, sigma_tr
  MaterialElastic<dim>::computeStressOnQuad(sigma_tr, grad_delta_u);
  sigma_tr += previous_sigma_tensor;

  // Compute deviatoric trial stress,  sigma_tr_dev
  sigma_tr_dev  = sigma_tr;
  sigma_tr_dev -= Matrix<Real>::eye(dim, sigma_tr.trace() / 3.0);

  // Compute effective deviatoric trial stress
  Real s = sigma_tr_dev.doubleDot(sigma_tr_dev);
  Real sigma_tr_dev_eff=std::sqrt(3./2. * s);

  // Compute deviatoric stress,  sigma_dev
  sigma_dev  = sigma_tensor;
  sigma_dev -= Matrix<Real>::eye(dim, sigma_tensor.trace() / 3.0);

  // Compute effective deviatoric stress
  s =  sigma_dev.doubleDot(sigma_dev);
  Real sigma_dev_eff = std::sqrt(3./2. * s);

  Real xr = 0.0;
  if(sigma_tr_dev_eff > sigma_dev_eff * Math::getTolerance())
    xr = sigma_dev_eff / sigma_tr_dev_eff;

  Real q = 1.5 * (1. / (1. +  3. * this->mu  / h) - xr);

  for (UInt m = 0; m < rows; ++m) {
    UInt i = VoigtHelper<dim>::vec[m][0];
    UInt j = VoigtHelper<dim>::vec[m][1];

    for (UInt n = 0; n < cols; ++n) {
      UInt k = VoigtHelper<dim>::vec[m][0];
      UInt l = VoigtHelper<dim>::vec[m][1];

      if (((sigma_tr_dev_eff-iso_hardening-sigmay) > 0) && (xr > 0)) {
	tangent(m,n) =
	  2. * this->mu * q * (sigma_tr_dev (i,j) / sigma_tr_dev_eff) * (sigma_tr_dev (k,l) / sigma_tr_dev_eff) +
	  (i==k) * (j==l) * 2. * this->mu * xr +
	  (i==j) * (k==l) * (this->kpa - 2./3. * this->mu * xr);
	if ((m == n) && (m>=dim))
	  tangent(m, n) = tangent(m, n) - this->mu * xr;
      } else {
	tangent(m,n) =
	  (i==k) * (j==l) * 2. * this->mu +
	  (i==j) * (k==l) * this->lambda;
	if ((m==n) && (m>=dim))
	  tangent(m,n) = tangent(m,n) - this->mu;
      }
      //correct tangent stiffness for shear component
    }
  }
}
