/**
 * @file   material_vreepeerlings_tmpl.hh
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
template<UInt spatial_dimension, template <UInt> class MatParent>
MaterialVreePeerlings<spatial_dimension, MatParent>::MaterialVreePeerlings(SolidMechanicsModel & model,
                                             const ID & id)  :
  Material(model, id),
  MaterialVreePeerlingsParent(model, id),
  Kapa("Kapa",id),
  strain_rate_vreepeerlings("strain-rate-vreepeerlings", id),
  Full_dam_value_strain("fulldam-valstrain", id),
  Full_dam_value_strain_rate("fulldam-valstrain-rate", id),
  Number_damage("number-damage", id),
  equi_strain("equi-strain", id),
  equi_strain_rate("equi-strain-rate", id)
 {
  AKANTU_DEBUG_IN();

  this->registerParam("Kapaoi"          , Kapaoi          , 0.0001, _pat_parsable);
  this->registerParam("Kapac"           , Kapac           , 0.0002, _pat_parsable);
  this->registerParam("Arate"           , Arate           , 0.    , _pat_parsable);
  this->registerParam("Brate"           , Brate           , 1.    , _pat_parsable);
  this->registerParam("Crate"           , Brate           , 1.    , _pat_parsable);
  this->registerParam("Kct"             , Kct             , 1.    , _pat_parsable);
  this->registerParam("Kapao_randomness", Kapao_randomness, 0.    , _pat_parsable);
  this->registerParam("Wf_vrplgs_ko_min", Wf_vrplgs_ko_min, 1.    , _pat_parsable);
  this->registerParam("Wf_vrplgs_lambda", Wf_vrplgs_lambda, 1.    , _pat_parsable);
  this->registerParam("Wf_vrplgs_mw"    , Wf_vrplgs_mw    , 1.    , _pat_parsable);

  this->initInternalArray(this->Kapa, 1);
  this->initInternalArray(this->equi_strain, 1);
  this->initInternalArray(this->equi_strain_rate, 1);
  this->initInternalArray(this->Full_dam_value_strain, 1);
  this->initInternalArray(this->Full_dam_value_strain_rate, 1);
  this->initInternalArray(this->Number_damage, 1);
  this->initInternalArray(this->strain_rate_vreepeerlings, spatial_dimension * spatial_dimension);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class MatParent>
void MaterialVreePeerlings<spatial_dimension, MatParent>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialVreePeerlingsParent::initMaterial();

  this->resizeInternalArray(this->Kapa);
  this->resizeInternalArray(this->equi_strain);
  this->resizeInternalArray(this->equi_strain_rate);
  this->resizeInternalArray(this->Full_dam_value_strain);
  this->resizeInternalArray(this->Full_dam_value_strain_rate);
  this->resizeInternalArray(this->Number_damage);
  this->resizeInternalArray(this->strain_rate_vreepeerlings);


  const Mesh & mesh = this->model->getFEM().getMesh();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension);

  for(; it != last_type; ++it) {
    Array <Real>::iterator<Real> kapa_it  = Kapa(*it).begin();
    Array <Real>::iterator<Real> kapa_end = Kapa(*it).end();

    for(; kapa_it != kapa_end; ++kapa_it) {
      Real limit = Kapac;
      Real rand = std::pow(std::abs(std::log(1 - drand48()*0.9999)), 1.0/Wf_vrplgs_mw );
      Real rand_part = Wf_vrplgs_lambda * rand * Kapao_randomness;
      Real rand_final = Wf_vrplgs_ko_min + rand_part;
       if (rand_final < limit)
      	{
	  *kapa_it = rand_final;
	  	}
	  else 
	  {
	  *kapa_it = limit;
	  }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class MatParent>
void MaterialVreePeerlings<spatial_dimension, MatParent>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MaterialVreePeerlingsParent::computeStress(el_type, ghost_type);

  Real * dam = this->damage(el_type, ghost_type).storage();
  Real * equi_straint = equi_strain(el_type, ghost_type).storage();
  Real * equi_straint_rate = equi_strain_rate(el_type, ghost_type).storage();
  Real * Kapaq = Kapa(el_type, ghost_type).storage();
  Real * FullDam_Valstrain = Full_dam_value_strain(el_type, ghost_type).storage();
  Real * FullDam_Valstrain_rate = Full_dam_value_strain_rate(el_type, ghost_type).storage();
  Real * Nb_damage = Number_damage(el_type, ghost_type).storage();
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
		      *FullDam_Valstrain,
		      *FullDam_Valstrain_rate,
		      *Nb_damage);
  ++dam;
  ++equi_straint;
  ++equi_straint_rate;
  ++Kapaq;
  ++strain_rate_vrplgs_it;
  ++FullDam_Valstrain;
  ++FullDam_Valstrain_rate;
  ++Nb_damage;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

