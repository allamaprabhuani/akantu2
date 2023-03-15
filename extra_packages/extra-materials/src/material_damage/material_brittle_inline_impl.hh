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
inline void MaterialBrittle<spatial_dimension>::computeStressOnQuad(
    Matrix<Real> & grad_u, Matrix<Real> & grad_v, Matrix<Real> & sigma,
    Real & dam, Real & sigma_equivalent, Real & fracture_stress) {

  MaterialElastic<spatial_dimension>::computeStressOnQuad(grad_u, sigma);

  Real equiv_strain_rate = 0.;
  Real volume_change_rate = grad_v.trace();
  if (spatial_dimension == 2) {
    equiv_strain_rate += 2. / 3. * pow(volume_change_rate / 3., 2.);
  }

  for (Int i = 0; i < spatial_dimension; ++i)
    for (Int j = 0; j < spatial_dimension; ++j)
      equiv_strain_rate += 2. / 3. *
                           pow(0.5 * (grad_v(i, j) + grad_v(j, i)) -
                                   (i == j) * volume_change_rate / 3.,
                               2.);

  equiv_strain_rate = sqrt(equiv_strain_rate);

  fracture_stress = S_0;
  if (equiv_strain_rate > E_0)
    fracture_stress = A;

  Vector<Real> principal_stress(spatial_dimension);
  sigma.eig(principal_stress);
  sigma_equivalent = principal_stress(0);
  for (Int i = 1; i < spatial_dimension; ++i)
    sigma_equivalent = std::max(sigma_equivalent, principal_stress(i));

  if (!this->is_non_local) {
    computeDamageAndStressOnQuad(sigma, dam, sigma_equivalent, fracture_stress);
  }
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
inline void MaterialBrittle<spatial_dimension>::computeDamageAndStressOnQuad(
    Matrix<Real> & sigma, Real & dam, Real & sigma_c, Real & fracture_stress) {
  if (sigma_c > fracture_stress)
    dam = 1.;

  dam = std::min(dam, 1.);

  sigma *= 1 - dam;
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
inline UInt MaterialBrittle<spatial_dimension>::getNbData(
    const Array<Element> & elements, const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = MaterialDamage<spatial_dimension>::getNbData(elements, tag);

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
inline void MaterialBrittle<spatial_dimension>::packData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  MaterialDamage<spatial_dimension>::packData(buffer, elements, tag);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
inline void
MaterialBrittle<spatial_dimension>::unpackData(CommunicationBuffer & buffer,
                                               const Array<Element> & elements,
                                               const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  MaterialDamage<spatial_dimension>::unpackData(buffer, elements, tag);

  AKANTU_DEBUG_OUT();
}
