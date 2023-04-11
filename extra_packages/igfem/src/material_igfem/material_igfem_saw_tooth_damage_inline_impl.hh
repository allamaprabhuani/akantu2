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
inline void
MaterialIGFEMSawToothDamage<spatial_dimension>::computeDamageAndStressOnQuad(
    Matrix<Real> & sigma, Real & dam) {
  sigma *= 1 - dam;
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
UInt MaterialIGFEMSawToothDamage<spatial_dimension>::updateDamage(
    UInt quad_index, const Real eq_stress, ElementType el_type,
    GhostType ghost_type) {
  AKANTU_DEBUG_ASSERT(prescribed_dam > 0.,
                      "Your prescribed damage must be greater than zero");

  Array<Real> & dam = this->damage(el_type, ghost_type);
  Real & dam_on_quad = dam(quad_index);

  /// check if damage occurs
  if (equivalent_stress(el_type, ghost_type)(quad_index) >=
      (1 - dam_tolerance) * norm_max_equivalent_stress) {
    /// damage the entire sub-element -> get the element index
    UInt el_index =
        quad_index / this->element_filter(el_type, ghost_type).getSize();
    UInt nb_quads = this->fem->getNbIntegrationPoints(el_type, ghost_type);
    UInt start_idx = el_index * nb_quads;
    Array<Idx> & sub_mat = this->sub_material(el_type, ghost_type);
    UInt damaged_quads = 0;
    if (dam_on_quad < dam_threshold) {
      for (Int q = 0; q < nb_quads; ++q, ++start_idx) {
        if (sub_mat(start_idx)) {
          dam(start_idx) += prescribed_dam;
          damaged_quads += 1;
        }
      }
    } else {
      for (Int q = 0; q < nb_quads; ++q, ++start_idx) {
        if (sub_mat(start_idx)) {
          dam(start_idx) += max_damage;
          damaged_quads += 1;
        }
      }
    }
    return damaged_quads;
  }

  return 0;
}

/* -------------------------------------------------------------------------- */
