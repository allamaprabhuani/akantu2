/**
 * @file   material_vreepeerlings_non_local_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date   Fri Jun 15 13:48:37 2012
 *
 * @brief  Specialization of the material class for the non-local Vree-Peerlings material
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

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeigthFunction, template <UInt> class MatParent>
MaterialVreePeerlingsNonLocal<spatial_dimension, WeigthFunction, MatParent>::MaterialVreePeerlingsNonLocal(SolidMechanicsModel & model,
												const ID & id)  :
  Material(model, id),
  MaterialVreePeerlingsNonLocalParent(model, id),
  equi_strain_non_local("equi-strain_non_local", id),
  equi_strain_rate_non_local("equi-strain-rate_non_local", id) {
  AKANTU_DEBUG_IN();

  this->is_non_local = true;

  this->initInternalArray(this->equi_strain_non_local, 1);
  this->initInternalArray(this->equi_strain_rate_non_local, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeigthFunction, template <UInt> class MatParent>
void MaterialVreePeerlingsNonLocal<spatial_dimension, WeigthFunction, MatParent>::initMaterial() {
  AKANTU_DEBUG_IN();

  this->resizeInternalArray(this->equi_strain_non_local);
  this->resizeInternalArray(this->equi_strain_rate_non_local);

  this->registerNonLocalVariable(this->equi_strain, this->equi_strain_non_local, 1);
  this->registerNonLocalVariable(this->equi_strain_rate, this->equi_strain_rate_non_local, 1);

  MaterialVreePeerlingsNonLocalParent::initMaterial();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

//template<UInt spatial_dimension, template <UInt> class WeigthFunction, template <UInt> class MatParent>
//void MaterialVreePeerlingsNonLocal<spatial_dimension, WeigthFunction, MatParent>::computeStress(ElementType el_type,
//										     GhostType ghost_type) {
//   AKANTU_DEBUG_IN();
//
//  Real * dam = this->damage(el_type, ghost_type).storage();
//  Real * equi_straint = equi_strain(el_type, ghost_type).storage();
//  Real * equi_straint_rate = equi_strain_rate(el_type, ghost_type).storage();
//  Real * Kapaq = this->Kapa(el_type, ghost_type).storage();
//  Real * crit_strain = this->critical_strain(el_type, ghost_type).storage();
//  Real * crit_strain_rate = this->critical_strain_rate(el_type, ghost_type).storage();
//  Real * rdr_damage = this->recorder_damage(el_type, ghost_type).storage();
//  Real  * nb_damage = this->number_damage(el_type, ghost_type).storage();
//  Real dt = this->model->getTimeStep();
//
//  Vector<UInt> & elem_filter = this->element_filter(el_type, ghost_type);
//  Vector<Real> & velocity = this->model->getVelocity();
//  Vector<Real> & strain_rate_vrplgs = this->strain_rate_vreepeerlings(el_type, ghost_type);
//
//
//  this->model->getFEM().gradientOnQuadraturePoints(velocity, strain_rate_vrplgs,
//						   spatial_dimension,
//						   el_type, ghost_type, &elem_filter);
//
//  Vector<Real>::iterator<types::RMatrix> strain_rate_vrplgs_it =
//    strain_rate_vrplgs.begin(spatial_dimension, spatial_dimension);
//
//
//  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
//
//  types::RMatrix & strain_rate = *strain_rate_vrplgs_it;
//
//
//
//  MaterialVreePeerlings<spatial_dimension>::computeStressOnQuad(grad_u, sigma,
//								*dam,
//								*equi_straint,
//								*equi_straint_rate,
//								*Kapaq,
//								dt,
//								strain_rate,
//								*crit_strain,
//								*crit_strain_rate,
//								*rdr_damage,
//								*nb_damage);
//  ++dam;
//  ++equi_straint;
//  ++equi_straint_rate;
//  ++Kapaq;
//  ++strain_rate_vrplgs_it;
//  ++crit_strain;
//  ++crit_strain_rate;
//  ++rdr_damage;
//  ++nb_damage;
//
//  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
//
//  AKANTU_DEBUG_OUT();
//}
//

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeigthFunction, template <UInt> class MatParent>
void MaterialVreePeerlingsNonLocal<spatial_dimension, WeigthFunction, MatParent>::computeNonLocalStress(ElementType el_type,
											     GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam   = this->damage(el_type, ghost_type).storage();
  Real * Kapaq = this->Kapa(el_type, ghost_type).storage();
  Real * equi_strain_nl = this->equi_strain_non_local(el_type, ghost_type).storage();
  Real * equi_strain_rate_nl = this->equi_strain_rate_non_local(el_type, ghost_type).storage();
  //Real * equi_strain_rate_nl = this->equi_strain_rate(el_type, ghost_type).storage();

  Real dt = this->model->getTimeStep();
  Real * FullDam_Valstrain = this->Full_dam_value_strain(el_type, ghost_type).storage();
  Real * FullDam_Valstrain_rate = this->Full_dam_value_strain_rate(el_type, ghost_type).storage();
  Real  * Nb_damage = this->Number_damage(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  this->computeDamageAndStressOnQuad(sigma, *dam, *equi_strain_nl, *equi_strain_rate_nl, *Kapaq, dt, *FullDam_Valstrain, *FullDam_Valstrain_rate, *Nb_damage);
  ++dam;
  ++equi_strain_nl;
  ++equi_strain_rate_nl;
  ++Kapaq;
  ++FullDam_Valstrain;
  ++FullDam_Valstrain_rate;
  ++Nb_damage;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

