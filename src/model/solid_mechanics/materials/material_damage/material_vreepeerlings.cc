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
  MaterialVreePeerlingsParent(model, id),
  Kapa("Kapa",id),
  strain_rate_vreepeerlings("strain-rate-vreepeerlings", id),
  critical_strain("critical-strain", id),
  critical_strain_rate("critical-strain-rate", id),
  recorder_damage("recorder-damage", id),
  number_damage("number-damage", id),
  equi_strain("equi-strain", id),
  equi_strain_rate("equi-strain-rate", id)
 {
  AKANTU_DEBUG_IN();

  this->registerParam("Kapa0i"          , Kapa0i          , 0.0001, _pat_parsable);
  this->registerParam("Kapa0"           , Kapa0           , 0.0001, _pat_parsable);
  this->registerParam("Alpha"           , Alpha           , 0.99  , _pat_parsable);
  this->registerParam("Beta"            , Beta            , 300.  , _pat_parsable);
  this->registerParam("Kct"             , Kct             , 1.    , _pat_parsable);
  this->registerParam("Kapa0_randomness", Kapa0_randomness, 0.    , _pat_parsable);

  // firststep = true;
  //countforinitialstep= 0;

  this->initInternalArray(this->Kapa, 1);
  this->initInternalArray(this->equi_strain, 1);
  this->initInternalArray(this->equi_strain_rate, 1);
  this->initInternalArray(this->critical_strain, 1);
  this->initInternalArray(this->critical_strain_rate, 1);
  this->initInternalArray(this->recorder_damage, 1);
  this->initInternalArray(this->number_damage, 1);
  this->initInternalArray(this->strain_rate_vreepeerlings, spatial_dimension * spatial_dimension);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialVreePeerlings<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialVreePeerlingsParent::initMaterial();

  this->resizeInternalArray(this->Kapa);
  this->resizeInternalArray(this->equi_strain);
  this->resizeInternalArray(this->equi_strain_rate);
  this->resizeInternalArray(this->critical_strain);
  this->resizeInternalArray(this->critical_strain_rate);
  this->resizeInternalArray(this->recorder_damage);
  this->resizeInternalArray(this->number_damage);
  this->resizeInternalArray(this->strain_rate_vreepeerlings);


  const Mesh & mesh = this->model->getFEM().getMesh();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension);

  Real beta = 1./3.;
  Real Lambda = 200e6 / this->E;

  for(; it != last_type; ++it) {
    Array <Real>::iterator<Real> kapa_it  = Kapa(*it).begin();
    Array <Real>::iterator<Real> kapa_end = Kapa(*it).end();

    for(; kapa_it != kapa_end; ++kapa_it) {
      // Real rand_part = (Kapa0i + drand48()*(Alpha - Kapa0i))*Kapa0_randomness;
      // Real rand_part = (2 * drand48()-1) * Kapa0_randomness * Kapa0i;
      Real rand = std::pow(std::abs(std::log(1 - drand48()*0.9999)), beta);
      Real rand_part = Lambda * rand * Kapa0_randomness;
      Real rand_final = Kapa0i + rand_part;
      if (rand_final < Alpha*.9)
	{
	  *kapa_it = rand_final;
	}
      else 
	{
	  *kapa_it = Alpha*.9;
	}
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialVreePeerlings<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MaterialVreePeerlingsParent::computeStress(el_type, ghost_type);

  Real * dam = this->damage(el_type, ghost_type).storage();
  Real * equi_straint = equi_strain(el_type, ghost_type).storage();
  Real * equi_straint_rate = equi_strain_rate(el_type, ghost_type).storage();
  Real * Kapaq = Kapa(el_type, ghost_type).storage();
  Real * crit_strain = critical_strain(el_type, ghost_type).storage();
  Real * crit_strain_rate = critical_strain_rate(el_type, ghost_type).storage();
  Real * rdr_damage = recorder_damage(el_type, ghost_type).storage();
  Real * nb_damage = number_damage(el_type, ghost_type).storage();
  Real dt = this->model->getTimeStep();

  Array<UInt> & elem_filter = this->element_filter(el_type, ghost_type);
  Array<Real> & velocity = this->model->getVelocity();
  Array<Real> & strain_rate_vrplgs = this->strain_rate_vreepeerlings(el_type, ghost_type);

  this->model->getFEM().gradientOnQuadraturePoints(velocity, strain_rate_vrplgs,
						   spatial_dimension,
						   el_type, ghost_type, elem_filter);

  Array<Real>::iterator< Matrix<Real> > strain_rate_vrplgs_it =
        strain_rate_vrplgs.begin(spatial_dimension, spatial_dimension);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  Matrix<Real> & strain_rate = *strain_rate_vrplgs_it;

  computeStressOnQuad(grad_u, sigma, 
		      *dam, 
		      *equi_straint,
		      *equi_straint_rate,
		      *Kapaq,
		      dt,
		      strain_rate,
		      *crit_strain,
		      *crit_strain_rate,
		      *rdr_damage,
		      *nb_damage);
  ++dam;
  ++equi_straint;
  ++equi_straint_rate;
  ++Kapaq;
  ++strain_rate_vrplgs_it;
  ++crit_strain;
  ++crit_strain_rate;
  ++rdr_damage;
  ++nb_damage;

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
//     const Array<Real> & kapa = Kapa(element.type, _not_ghost);
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
//     Array<Real> & kapa = Kapa(element.type, _not_ghost);
//     for(UInt q = 0; q < nb_quad; ++q)
//       buffer >> kapa(element.element * nb_quad + q);
//   }

//   MaterialDamage<spatial_dimension>::packData(buffer, element, tag);
//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialVreePeerlings);



__END_AKANTU__
