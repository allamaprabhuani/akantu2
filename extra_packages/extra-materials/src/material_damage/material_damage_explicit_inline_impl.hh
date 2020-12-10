/**
 * @file   material_damage_explicit_inline_impl.hh
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @date   Tue Feb 12 2019
 *
 * @brief  Inline implementation of material damage explicit
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
#include "communicator.hh"
#include "material_damage_explicit.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim, template <UInt> class ElasticParent>
MaterialDamageExplicit<dim, ElasticParent>::MaterialDamageExplicit(
    SolidMechanicsModel & model, const ID & id)
    : parent(model, id), Sc("Sc", *this), Gf(0.), crack_band_width(0.),
      max_damage(1.), Et("Et", *this) {

  AKANTU_DEBUG_IN();
  this->registerParam("Sc", Sc, _pat_parsable, "critical stress threshold");
  this->registerParam("Gf", Gf, _pat_parsable | _pat_modifiable,
                      "fracture energy");
  this->registerParam("crack_band_width", crack_band_width,
                      _pat_parsable | _pat_modifiable, "crack_band_width");
  this->registerParam("max_damage", max_damage, _pat_parsable | _pat_modifiable,
                      "maximum possible damage");

  this->Sc.initialize(1);
  this->Et.initialize(1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim, template <UInt> class ElasticParent>
void MaterialDamageExplicit<dim, ElasticParent>::initMaterial() {
  AKANTU_DEBUG_IN();

  parent::initMaterial();

  GhostType ghost_type = _not_ghost;

  for (auto && el_type : this->element_filter.elementTypes(dim, ghost_type)) {
    for (auto && data :
         zip(this->Sc(el_type, ghost_type), this->Et(el_type, ghost_type))) {
      auto & strength = std::get<0>(data);
      auto & soft_slope = std::get<1>(data);
      soft_slope = strength / (strength / this->E - 2 * this->Gf / strength /
                                                        this->crack_band_width);
      AKANTU_DEBUG_ASSERT(soft_slope < 0,
                          "The softening slope is positive or 0");
    }
  }
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt dim, template <UInt> class ElasticParent>
void MaterialDamageExplicit<dim, ElasticParent>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  parent::computeStress(el_type, ghost_type);

  auto dam = this->damage(el_type, ghost_type).begin();
  auto Et = this->Et(el_type, ghost_type).begin();
  auto Sc = this->Sc(el_type, ghost_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  computeDamageAndStressOnQuad(grad_u, sigma, *dam, *Et, *Sc);

  ++dam;
  ++Et;
  ++Sc;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  /*  computeNormalizedEquivalentStress(el_type, ghost_type);
    updateSofteningStrain(el_type, ghost_type);
    this->norm_max_equivalent_stress = 0;
    findMaxAndAvNormalizedEquivalentStress(el_type, ghost_type);
  */
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
//  template <UInt dim, template <UInt> class ElasticParent>
// void MaterialDamageExplicit<dim, ElasticParent>::
//     computeNormalizedEquivalentStress(ElementType el_type,
//                                       GhostType ghost_type) {
//   AKANTU_DEBUG_IN();

//   /// Vector to store eigenvalues of current stress tensor
//   Vector<Real> eigenvalues(dim);

//   for (auto && data :
//        zip(make_view(this->stress(el_type, ghost_type), dim, dim),
//            this->Sc(el_type, ghost_type),
//            this->equivalent_stress(el_type, ghost_type),
//            this->eig_pos(el_type, ghost_type),
//            make_view(this->rot_matrix(el_type, ghost_type), dim, dim))) {

//     const auto & sigma = std::get<0>(data);
//     const auto & sigma_crit = std::get<1>(data);
//     auto & sigma_equ = std::get<2>(data);
//     auto & pos = std::get<3>(data);
//     auto & rotation_matrix = std::get<4>(data);

//     /// compute eigenvalues
//     sigma.eig(eigenvalues, rotation_matrix, false);

//     /// find max eigenvalue and normalize by tensile strength
//     sigma_equ = *(std::max_element(eigenvalues.storage(),
//                                    eigenvalues.storage() + dim)) /
//                 sigma_crit;

//     /// normalize each column of the rotation matrix by the length of
//     /// corresponding eigen vector
//     for (auto && c : arange(rotation_matrix.cols())) {
//       Vector<Real> vect = rotation_matrix(c);
//       vect /= vect.norm();
//     }

//     pos = std::distance(
//         eigenvalues.storage(),
//         std::max_element(eigenvalues.storage(), eigenvalues.storage() +
//         dim));
//   }

//   AKANTU_DEBUG_OUT();
// }

// /* --------------------------------------------------------------------------
// */ template <UInt dim, template <UInt> class ElasticParent> void
// MaterialDamageExplicit<dim, ElasticParent>::updateSofteningStrain(
//     ElementType el_type, GhostType ghost_type) {
//   AKANTU_DEBUG_IN();
//   if (this->getID().find("damage_viscoelastic") != std::string::npos) {
//     for (auto && data :
//          zip(make_view(this->gradu(el_type, ghost_type), dim, dim),
//              this->equivalent_stress(el_type, ghost_type),
//              this->eig_pos(el_type, ghost_type),
//              make_view(this->rot_matrix(el_type, ghost_type), dim, dim),
//              this->eps_0(el_type, ghost_type),
//              this->damage(el_type, ghost_type))) {

//       const auto & strain = std::get<0>(data);
//       const auto & sigma_equ = std::get<1>(data);
//       const auto & pos = std::get<2>(data);
//       const auto & rotation_matrix = std::get<3>(data);
//       auto & eps_soft = std::get<4>(data);
//       const auto & dam = std::get<5>(data);

//       if (!dam) {
//         /// check if softening starts
//         if (sigma_equ >= (1 - this->dam_tolerance)) {
//           Matrix<Real> eig_strain(dim, dim);
//           Matrix<Real> temp(eig_strain);
//           temp.mul<false, false>(rotation_matrix, strain);
//           eig_strain.mul<false, true>(temp, rotation_matrix);
//           eps_soft = eig_strain(pos, pos);
//         }
//       }
//       AKANTU_DEBUG_OUT();
//     }
//   }
// }

// /* -------------------------------------------------------------------------
// */ template <UInt dim, template <UInt> class ElasticParent> UInt
// MaterialDamageExplicit<dim, ElasticParent>::updateDamage() {
//   UInt nb_damaged_elements = 0;

//   if (this->norm_max_equivalent_stress >= 1.) {

//     AKANTU_DEBUG_IN();

//     /// update the damage only on non-ghosts elements! Doesn't make sense to
//     /// update on ghost.
//     GhostType ghost_type = _not_ghost;

//     for (auto && el_type : this->element_filter.elementTypes(dim,
//     ghost_type)) {

//       /// loop over all the quads of the given element type
//       for (auto && data :
//            zip(this->equivalent_stress(el_type, ghost_type),
//                this->damage(el_type, ghost_type),
//                this->eps_0(el_type, ghost_type), this->Sc(el_type,
//                ghost_type), make_view(this->gradu(el_type, ghost_type), dim,
//                dim), this->eig_pos(el_type, ghost_type),
//                make_view(this->rot_matrix(el_type, ghost_type), dim, dim))) {
//         auto & eq_stress = std::get<0>(data);
//         auto & damage = std::get<1>(data);
//         auto & eps_0 = std::get<2>(data);
//         auto & strength = std::get<3>(data);
//         const auto & strain = std::get<4>(data);
//         const auto & pos = std::get<5>(data);
//         const auto & rotation_matrix = std::get<6>(data);

//         /// check if damage occurs
//         if (eq_stress >=
//             (1 - this->dam_tolerance) * this->norm_max_equivalent_stress) {

//           /// principal strain corresponding to the largest principal stress
//           Matrix<Real> eig_strain(dim, dim);
//           Matrix<Real> temp(eig_strain);
//           temp.mul<false, false>(rotation_matrix, strain);
//           eig_strain.mul<false, true>(temp, rotation_matrix);
//           Real eps = eig_strain(pos, pos);

//           /// update the damage on this quad based on the strain
//           damage = std::max(damage, (eps - eps_0) / (2 * this->Gf / strength
//           /
//                                                      this->crack_band_width));
//           damage = std::min(damage, this->max_damage);
//           strength *= (1 - damage);
//           if (!damage)
//             nb_damaged_elements += 1;
//         }
//       }
//     }
//   }

//   AKANTU_DEBUG_OUT();
//   return nb_damaged_elements;
// }

/* -------------------------------------------------------------------------- */
template <UInt dim, template <UInt> class ElasticParent>
inline void
MaterialDamageExplicit<dim, ElasticParent>::computeDamageAndStressOnQuad(
    Matrix<Real> & grad_u, Matrix<Real> & sigma, Real & dam, const Real & Et,
    const Real & Sc) {
  Real w = 0;
  for (UInt i = 0; i < dim; ++i) {
    for (UInt j = 0; j < dim; ++j) {
      w += sigma(i, j) * (grad_u(i, j) + grad_u(j, i)) / 2.;
    }
  }
  w *= 0.5;
  //w *= (1 - dam);

  Real wy = Sc * Sc / (2 * this->E);
  Real gamma = -Et / this->E;
  Real k = wy * pow(((1 + gamma) / (1 + gamma - dam)), 2);

  Real F = w - k;

  if (F > 0) {
    dam = (1 + gamma) * (1 - std::sqrt(wy / w));
    dam = std::min(dam, this->max_damage);
    dam = std::max(0., dam);
  }
  sigma *= 1 - dam;
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
