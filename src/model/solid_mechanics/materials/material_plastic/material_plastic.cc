/**
 * @file   material_plastic.cc
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Fri Apr  4 13:37:16 2014
 *
 * @brief
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
#include "material_plastic.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialPlastic<spatial_dimension>::MaterialPlastic(SolidMechanicsModel & model, const ID & id) :
  Material(model, id),
  MaterialElastic<spatial_dimension>(model, id),
  iso_hardening("iso_hardening", *this),
  inelastic_strain("inelastic_strain", *this),
  plastic_energy("plastic_energy", *this),
  d_plastic_energy("d_plastic_energy", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam("h", h, 0., _pat_parsable | _pat_modifiable, "Hardening  modulus");
  this->registerParam("sigma_y", sigma_y, 0., _pat_parsable | _pat_modifiable, "Yield stress");

  this->iso_hardening.initialize(1);
  this->iso_hardening.initializeHistory();

  this->plastic_energy.initialize(1);
  this->d_plastic_energy.initialize(1);

  this->finite_deformation            = false;
  this->use_previous_stress           = true;
  this->use_previous_strain           = true;
  this->use_previous_stress_thermal   = true;


  this->inelastic_strain.initialize(spatial_dimension * spatial_dimension);
  this->inelastic_strain.initializeHistory();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialPlastic<spatial_dimension>::getEnergy(std::string type) {
  AKANTU_DEBUG_IN();

  if (type == "plastic") return getPlasticEnergy();
  else return MaterialElastic<spatial_dimension>::getEnergy(type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialPlastic<spatial_dimension>::getPlasticEnergy() {
  AKANTU_DEBUG_IN();

  Real penergy = 0.;

  const Mesh & mesh = this->model->getFEM().getMesh();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _not_ghost);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, _not_ghost);

  for(; it != end; ++it) {
    penergy += this->model->getFEM().integrate(plastic_energy(*it, _not_ghost),
                                               *it, _not_ghost,
                                               this->element_filter(*it, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return penergy;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialPlastic<spatial_dimension>::computePotentialEnergy(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if(ghost_type != _not_ghost) return;

  Array<Real>::scalar_iterator epot = this->potential_energy(el_type, ghost_type).begin();

   Array<Real>::const_iterator< Matrix<Real> > inelastic_strain_it
     = this->inelastic_strain(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  Matrix<Real> elastic_strain(spatial_dimension, spatial_dimension);
  elastic_strain.copy(grad_u);
  elastic_strain -= *inelastic_strain_it;

  MaterialElastic<spatial_dimension>::computePotentialEnergyOnQuad(elastic_strain, sigma, *epot);

  ++epot;
  ++inelastic_strain_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialPlastic<spatial_dimension>::updateEnergies(ElementType el_type,
                                                        GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MaterialElastic<spatial_dimension>::updateEnergies(el_type, ghost_type);

  Array<Real>::iterator<> pe_it =
    this->plastic_energy(el_type, ghost_type).begin();

  Array<Real>::iterator<> wp_it =
    this->d_plastic_energy(el_type, ghost_type).begin();

  Array<Real>::iterator< Matrix<Real> > inelastic_strain_it =
    this->inelastic_strain(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

  Array<Real>::iterator< Matrix<Real> > previous_inelastic_strain_it =
    this->inelastic_strain.previous(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

  Array<Real>::matrix_iterator previous_sigma =
    this->stress.previous(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);


  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  Matrix<Real> delta_strain_it(*inelastic_strain_it);
  delta_strain_it -= *previous_inelastic_strain_it;

  Matrix<Real> sigma_h(sigma);
  sigma_h += *previous_sigma;

  *wp_it = .5 * sigma_h.doubleDot(delta_strain_it);

  *pe_it += *wp_it;

  ++pe_it;
  ++wp_it;
  ++inelastic_strain_it;
  ++previous_inelastic_strain_it;
  ++previous_sigma;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialPlastic);

__END_AKANTU__
