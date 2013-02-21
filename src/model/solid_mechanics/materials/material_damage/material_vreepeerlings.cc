/**
 * @file   material_vreepeerlings.cc
 *
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date   Fri Feb 24 14:27:15 2012
 *
 * @brief  Specialization of the material class for the VreePeerlings material
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
#include "material_vreepeerlings.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialVreePeerlings<spatial_dimension>::MaterialVreePeerlings(SolidMechanicsModel & model,
					     const ID & id)  :
  Material(model, id),
  MaterialDamage<spatial_dimension>(model, id),
  Kapa("Kapa",id),
  strain_rate_vreepeerlings("strain-rate-vreepeerlings", id),
  critical_strain("critical-strain", id)
 {
  AKANTU_DEBUG_IN();

  this->registerParam("Kapa0i"          , Kapa0i          , 0.0001, _pat_parsable);
  this->registerParam("Kapa0"           , Kapa0           , 0.0001, _pat_parsable);
  this->registerParam("Alpha"           , Alpha           , 0.99  , _pat_parsable);
  this->registerParam("Beta"            , Beta            , 300.  , _pat_parsable);
  this->registerParam("Kct"             , Kct             , 1.    , _pat_parsable);
  this->registerParam("Kapa0_randomness", Kapa0_randomness, 0.    , _pat_parsable);

  firststep = true;
  countforinitialstep= 0;

  this->initInternalVector(this->Kapa, 1);
  this->initInternalVector(this->critical_strain, 1);
  this->initInternalVector(this->strain_rate_vreepeerlings, spatial_dimension * spatial_dimension);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialVreePeerlings<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialDamage<spatial_dimension>::initMaterial();

  this->resizeInternalVector(this->Kapa);
  this->resizeInternalVector(this->critical_strain);
  this->resizeInternalVector(this->strain_rate_vreepeerlings);

  const Mesh & mesh = this->model->getFEM().getMesh();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension);

  for(; it != last_type; ++it) {
    Vector <Real>::iterator<Real> kapa_it  = Kapa(*it).begin();
    Vector <Real>::iterator<Real> kapa_end = Kapa(*it).end();

    for(; kapa_it != kapa_end; ++kapa_it) {
      Real rand_part = (2 * drand48()-1) * Kapa0_randomness * Kapa0i;
      *kapa_it = Kapa0i + rand_part;
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialVreePeerlings<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam = this->damage(el_type, ghost_type).storage();
  Real * Kapaq = Kapa(el_type, ghost_type).storage();
  Real * crit_strain = critical_strain(el_type, ghost_type).storage();
  Real dt = this->model->getTimeStep();

  Vector<UInt> & elem_filter = this->element_filter(el_type, ghost_type);
  Vector<Real> & velocity = this->model->getVelocity();
  Vector<Real> & strain_rate_vrplgs = this->strain_rate_vreepeerlings(el_type, ghost_type); 

  this->model->getFEM().gradientOnQuadraturePoints(velocity, strain_rate_vrplgs,
						   spatial_dimension,
						   el_type, ghost_type, &elem_filter);

  Vector<Real>::iterator<types::RMatrix> strain_rate_vrplgs_it =
        strain_rate_vrplgs.begin(spatial_dimension, spatial_dimension);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  Real Equistrain;
  Real Equistrain_rate;
  types::RMatrix & strain_rate = *strain_rate_vrplgs_it;

  computeStressOnQuad(grad_u, sigma, *dam, Equistrain, Equistrain_rate, *Kapaq, dt, strain_rate, *crit_strain);
  ++dam;
  ++Kapaq;
  ++strain_rate_vrplgs_it;
  ++crit_strain;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
//UInt MaterialVreePeerlings<spatial_dimension>::getNbDataToPack(const Element & element,
//					    SynchronizationTag tag) const {
//  AKANTU_DEBUG_IN();
//   UInt size = 0;
//   if(tag == _gst_smm_init_mat) {
//     UInt nb_quad = this->model->getFEM().getNbQuadraturePoints(element.type);
//     size += sizeof(Real) + nb_quad;
//   }

//   size += MaterialDamage<spatial_dimension>::getNbDataToPack(element, tag);

//   AKANTU_DEBUG_OUT();
//   return size;
// }

// /* -------------------------------------------------------------------------- */
// template<UInt spatial_dimension>
// UInt MaterialVreePeerlings<spatial_dimension>::getNbDataToUnpack(const Element & element,
// 					      SynchronizationTag tag) const {
//   AKANTU_DEBUG_IN();

//   UInt size = 0;
//   if(tag == _gst_smm_init_mat) {
//     UInt nb_quad = this->model->getFEM().getNbQuadraturePoints(element.type);
//     size += sizeof(Real) + nb_quad;
//   }

//   size += MaterialDamage<spatial_dimension>::getNbDataToPack(element, tag);

//   AKANTU_DEBUG_OUT();
//   return size;
// }

// /* -------------------------------------------------------------------------- */
// template<UInt spatial_dimension>
// void MaterialVreePeerlings<spatial_dimension>::packData(CommunicationBuffer & buffer,
// 				     const Element & element,
// 				     SynchronizationTag tag) const {
//   AKANTU_DEBUG_IN();

//   if(tag == _gst_smm_init_mat){
//     UInt nb_quad = this->model->getFEM().getNbQuadraturePoints(element.type);
//     const Vector<Real> & kapa = Kapa(element.type, _not_ghost);
//     for(UInt q = 0; q < nb_quad; ++q)
//       buffer << kapa(element.element * nb_quad + q);
//   }

//   MaterialDamage<spatial_dimension>::packData(buffer, element, tag);
//   AKANTU_DEBUG_OUT();
// }

// /* -------------------------------------------------------------------------- */
// template<UInt spatial_dimension>
// void MaterialVreePeerlings<spatial_dimension>::unpackData(CommunicationBuffer & buffer,
// 				       const Element & element,
// 				       SynchronizationTag tag) {
//   AKANTU_DEBUG_IN();

//   if(tag == _gst_smm_init_mat) {
//     UInt nb_quad = this->model->getFEM().getNbQuadraturePoints(element.type);
//     Vector<Real> & kapa = Kapa(element.type, _not_ghost);
//     for(UInt q = 0; q < nb_quad; ++q)
//       buffer >> kapa(element.element * nb_quad + q);
//   }

//   MaterialDamage<spatial_dimension>::packData(buffer, element, tag);
//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
INSTANSIATE_MATERIAL(MaterialVreePeerlings);



__END_AKANTU__
