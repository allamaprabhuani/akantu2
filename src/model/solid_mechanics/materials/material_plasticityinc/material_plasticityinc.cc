/**
 * @file   material_plasticity_inc.cc
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date   Tue Jul 09 18:15:37 20130
 *
 * @brief  Specialization of the material class for isotropic finite deformation linear hardening plasticity 
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
#include "material_plasticityinc.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt dim>
MaterialPlasticityinc<dim>::MaterialPlasticityinc(SolidMechanicsModel & model, const ID & id) :
  Material(model, id),
  MaterialElastic<dim>(model, id),
  iso_hardening("iso_hardening", *this),
  plastic_energy("plastic_energy", *this),
  d_plastic_energy("d_plastic_energy", *this),
  d_inelas_strain("d_inelas_strain", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam("h", h, 0., _pat_parsable | _pat_modifiable, "Hardening  modulus");
  this->registerParam("sigmay", sigmay, 0., _pat_parsable | _pat_modifiable, "Yield stress");

  this->iso_hardening.initialize(1);
  this->plastic_energy.initialize(1);
  this->d_plastic_energy.initialize(1);
  this->d_inelas_strain.initialize(dim * dim);

  this->inelastic_deformation       = true;
  this->finite_deformation          = false;
  this->use_previous_stress         = true;
  this->use_previous_stress_thermal = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialPlasticityinc<dim>::computeAllStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Material::computeAllStresses(ghost_type);

  Mesh::type_iterator it = this->model->getFEM().getMesh().firstType(dim, ghost_type);
  Mesh::type_iterator last_type = this->model->getFEM().getMesh().lastType(dim, ghost_type);

  for(; it != last_type; ++it) {
    computePlasticEnergyIncrement(*it, ghost_type);
    updatePlasticEnergy(*it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialPlasticityinc<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MaterialThermal<spatial_dimension>::computeStress(el_type, ghost_type);

  Array<Real>::iterator<> sigma_th_cur_it = 
    this->sigma_th_cur(el_type, ghost_type).begin();

  Array<Real>::iterator<> sigma_th_prev_it = 
    this->sigma_th_prev(el_type, ghost_type).begin();

  Array<Real>::iterator< Matrix<Real> > d_strain_it =
    this->delta_strain(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

  Array<Real>::iterator< Matrix<Real> > inelas_strain_it =
    this->inelas_strain(el_type, ghost_type).begin(spatial_dimension,spatial_dimension);

  Array<Real>::iterator< Matrix<Real> > d_inelas_strain_it =
    this->d_inelas_strain(el_type, ghost_type).begin(spatial_dimension,spatial_dimension);

  Array<Real>::iterator<> iso_hardening_it =
    this->iso_hardening(el_type, ghost_type).begin();


  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  Matrix<Real> & grad_delta_u = *d_strain_it;
  Matrix<Real> & inelas_strain_tensor = *inelas_strain_it;
  Matrix<Real> & d_inelas_strain_tensor = *d_inelas_strain_it;

  computeStressOnQuad(grad_u,
                      grad_delta_u,
                      sigma,
                      inelas_strain_tensor,
                      d_inelas_strain_tensor,
                      *iso_hardening_it,
                      *sigma_th_cur_it,
                      *sigma_th_prev_it);
  ++d_strain_it;
  ++inelas_strain_it;
  ++d_inelas_strain_it;
  ++iso_hardening_it;
  ++sigma_th_cur_it;
  ++sigma_th_prev_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialPlasticityinc<spatial_dimension>::getEnergy(std::string type) {
  AKANTU_DEBUG_IN();

  if (type == "plastic") return getPlasticEnergy();
  else return MaterialElastic<spatial_dimension>::getEnergy(type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialPlasticityinc<spatial_dimension>::getPlasticEnergy() {
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
void MaterialPlasticityinc<spatial_dimension>::computePotentialEnergy(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if(ghost_type != _not_ghost) return;
  Real * epot = this->potential_energy(el_type, ghost_type).storage();

   Array<Real>::const_iterator< Matrix<Real> > inelas_strain_it = this->inelas_strain(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  computePotentialEnergyOnQuad(grad_u, sigma, *inelas_strain_it, *epot);
  ++epot;
  ++inelas_strain_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialPlasticityinc<spatial_dimension>::computeTangentModuli(__attribute__((unused)) const ElementType & el_type,
								    Array<Real> & tangent_matrix,
								    __attribute__((unused)) GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Int dim=spatial_dimension;

  Array<Real>::iterator< Matrix<Real> > previous_sigma_it =
    this->previous_stress(el_type, ghost_type).begin(dim, dim);

  Real * iso_hardening= this->iso_hardening(el_type, ghost_type).storage();


  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_FINITE_DEFORMATION_BEGIN(tangent_matrix);

  //   Matrix<Real> & grad_delta_u = *d_strain_it;
  //   Matrix<Real> & sigma_tensor = *sigma_it;
  Matrix<Real> & previous_sigma_tensor = *previous_sigma_it;

  computeTangentModuliOnQuad(tangent, grad_delta_u, sigma_tensor, previous_sigma_tensor, *iso_hardening);
  //   ++green_it; ++sigma_it;
  ++previous_sigma_it;
  ++iso_hardening;
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_FINITE_DEFORMATION_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialPlasticityinc<spatial_dimension>::updatePlasticEnergy(ElementType el_type,
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<Real>::iterator<> pe_it =
    this->plastic_energy(el_type, ghost_type).begin();
  Array<Real>::iterator<> wp_it =
    this->d_plastic_energy(el_type, ghost_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  *pe_it += *wp_it;

  ++pe_it;
  ++wp_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialPlasticityinc<spatial_dimension>::computePlasticEnergyIncrement(ElementType el_type,
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<Real>::iterator<> wp_it =
    this->d_plastic_energy(el_type, ghost_type).begin();
  Array<Real>::iterator< Matrix<Real> > d_strain_it =
    this->d_inelas_strain(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  *wp_it = sigma.doubleDot(*d_strain_it);
  ++wp_it;
  ++d_strain_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialPlasticityinc);

__END_AKANTU__
