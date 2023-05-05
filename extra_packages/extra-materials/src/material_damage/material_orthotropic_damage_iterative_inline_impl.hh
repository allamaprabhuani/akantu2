/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
inline void MaterialOrthotropicDamageIterative<spatial_dimension>::
    computeDamageAndStressOnQuad(Matrix<Real> & sigma,
                                 Matrix<Real> & one_minus_D,
                                 Matrix<Real> & sqrt_one_minus_D,
                                 Matrix<Real> & damage,
                                 Matrix<Real> & first_term,
                                 Matrix<Real> & third_term) {

  // Real dmax = *(std::max_element(damage.data(), damage.data() +
  // spatial_dimension*spatial_dimension) );
  Real eta_effective = 0;

  // if ( (1 - dmax*dmax)  < (1 - this->eta / spatial_dimension *
  // damage.trace()) ) {

  //   eta_effective = this->spatial_dimension * dmax * dmax / damage.trace();

  // }
  // else
  eta_effective = this->eta;

  /// hydrostatic sensitivity parameter
  // Real eta = 3.;

  /// Definition of Cauchy stress based on second order damage tensor:
  /// "Anisotropic damage modelling of biaxial behaviour and rupture
  /// of concrete strucutres", Ragueneau et al., 2008, Eq. 7
  first_term.mul<false, false>(sqrt_one_minus_D, sigma);
  first_term *= sqrt_one_minus_D;

  Real second_term = 0;
  for (Int i = 0; i < this->spatial_dimension; ++i) {
    for (Int j = 0; j < this->spatial_dimension; ++j)
      second_term += sigma(i, j) * one_minus_D(i, j);
  }

  second_term /= (this->spatial_dimension - damage.trace());

  // for (Int i = 0; i < this->spatial_dimension; ++i) {
  //   for (Int j = 0; j < this->spatial_dimension; ++j)
  //     one_minus_D(i,j) *= second_term;
  // }
  one_minus_D *= second_term;

  third_term.eye(
      1. / this->spatial_dimension * sigma.trace() *
      (1 - std::min(eta_effective / (this->spatial_dimension) * damage.trace(),
                    this->max_damage)));

  sigma.copy(first_term);
  sigma -= one_minus_D;
  sigma += third_term;
}
