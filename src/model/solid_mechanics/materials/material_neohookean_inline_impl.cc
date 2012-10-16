/**
 * @file   material_neohookean_inline_impl.cc
 * @author Marion Chambart <marion.chambart@epfl.ch>
 * @date   Tue Jul 27 11:57:43 2010
 *
 * @brief  Implementation of the inline functions of the material neohookean
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

/* -------------------------------------------------------------------------- */
#include "aka_math.hh"

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void
MaterialNeohookean<spatial_dimension>::computeStressOnQuad(types::RMatrix & grad_u,
							   types::RMatrix & sigma) {
  types::RMatrix F(3, 3);
  this->template gradUToF<spatial_dimension>(grad_u, F);

  ///First compute the Left Cauchy-Green deformation tensor : C= F^tF.
  types::RMatrix C(3, 3);
  this->rightCauchy(F, C);

  ///Compute determinant of C
  Real detC = Math::det3(C.storage());
  Real defvol = 0.5 * log(detC);

  Real p = lambda * defvol;

  types::RMatrix S(3, 3);
  Math::inv3(C.storage(), S.storage());

  S *= p - mu;

  for (UInt i = 0; i < 3; ++i) S(i,i) = S(i,i) + mu;

  types::RMatrix sigma_tmp(3, 3);
  sigma_tmp.mul<false, false>(F, S);

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      sigma(i, j) = sigma_tmp(i, j);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void
MaterialNeohookean<spatial_dimension>::computeTangentModuliOnQuad(types::RMatrix & grad_u,
								  types::RMatrix & tangent) {
  UInt n = tangent.cols();
  types::RMatrix F(3, 3);
  this->template gradUToF<spatial_dimension>(grad_u, F);
  Real J = Math::det3(F.storage());
  Real Miiii = 2*mu + lambda;
  Real Miijj = lambda*J*J;
  Real Mijij = mu - 0.5*lambda*(J*J - 1);

  tangent(0, 0) = Miiii;

  // test of dimension should by optimized out by the compiler due to the template
  if(spatial_dimension >= 2) {
    tangent(1, 1) = Miiii;
    tangent(0, 1) = Miijj;
    tangent(1, 0) = Miijj;

    tangent(n - 1, n - 1) = Mijij;
  }

  if(spatial_dimension == 3) {
    tangent(2, 2) = Miiii;
    tangent(0, 2) = Miijj;
    tangent(1, 2) = Miijj;
    tangent(2, 0) = Miijj;
    tangent(2, 1) = Miijj;

    tangent(3, 3) = Mijij;
    tangent(4, 4) = Mijij;
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void
MaterialNeohookean<spatial_dimension>::computePotentialEnergyOnQuad(types::RMatrix & grad_u,
								    Real & epot) {

  types::RMatrix F(3, 3);
  types::RMatrix C(3, 3);

  Material::gradUToF<spatial_dimension>(grad_u, F);
  this->rightCauchy(F, C);

  Real detC = Math::det3(C.storage());

  Real defvol = 0.5*log(detC);
  Real p = lambda * defvol;
  Real traceC = C.trace();

  /// potential energy
  epot = (0.5*p - mu)*defvol + 0.5*mu*(traceC - 3.);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline Real MaterialNeohookean<spatial_dimension>::getStableTimeStep(Real h,
								     const Element & element) {
  return (h/celerity(element));
}
