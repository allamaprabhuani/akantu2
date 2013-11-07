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
  MaterialThermal<dim>(model, id),
  iso_hardening("iso_hardening", id) {
  AKANTU_DEBUG_IN();

  this->registerParam("h", h, 0., _pat_parsable | _pat_modifiable, "Hardening  modulus");
  this->registerParam("sigmay", sigmay, 0., _pat_parsable | _pat_modifiable, "Yield stress");
  this->registerParam("Plane_Stress", plane_stress, false, _pat_parsmod, "Is plane stress"); /// @todo Plane_Stress should not be possible to be modified after initMaterial (but before)
  this->registerParam("lambda", lambda, _pat_readable, "First Lamé coefficient");
  this->registerParam("mu", mu, _pat_readable, "Second Lamé coefficient");
  this->registerParam("kapa", kpa, _pat_readable, "Bulk coefficient");

  this->initInternalArray(this->iso_hardening, 1);

  this->inelastic_deformation=true;
  this->finite_deformation=false;
  this->use_previous_stress=true;
  this->use_previous_stress_thermal=true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialPlasticityinc<dim>::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();
  MaterialThermal<dim>::initMaterial();
  this->resizeInternalArray(this->iso_hardening);
  if (dim == 1) this->nu = 0.;
  updateInternalParameters();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialPlasticityinc<dim>::updateInternalParameters() {
  MaterialThermal<dim>::updateInternalParameters();
  lambda = this->nu * this->E / ((1 + this->nu) * (1 - 2 * this->nu));
  mu = this->E / (2 * (1 + this->nu));

  if (plane_stress) lambda = this->nu * this->E / ((1 + this->nu)*(1 - this->nu));

  kpa = lambda + 2. / 3. * mu;
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialPlasticityinc<dim>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();


  //Array<UInt> & elem_filter = element_filter(el_type, ghost_type);
  //UInt nb_element = elem_filter.getSize();

  MaterialThermal<dim>::computeStress(el_type, ghost_type);

  Array<Real>::iterator<> sigma_th_cur_it = 
    this->sigma_th_cur(el_type, ghost_type).begin();

  Array<Real>::iterator<> sigma_th_prev_it = 
    this->sigma_th_prev(el_type, ghost_type).begin();

  Array<Real>::iterator< Matrix<Real> > sigma_it =
    this->stress(el_type, ghost_type).begin(dim,dim);

  Array<Real>::iterator< Matrix<Real> > strain_it =
    this->strain(el_type, ghost_type).begin(dim,dim);

  Array<Real>::iterator< Matrix<Real> > strain_end =
    this->strain(el_type, ghost_type).end(dim,dim);

  Array<Real>::iterator< Matrix<Real> > d_strain_it =
    this->delta_strain(el_type, ghost_type).begin(dim, dim);

  Array<Real>::iterator< Matrix<Real> > inelas_strain_it =
    this->inelas_strain(el_type, ghost_type).begin(dim,dim);

  Real * iso_hardening= this->iso_hardening(el_type, ghost_type).storage();

  //Array<Real>::iterator< Matrix<Real> > previous_stress_it =
  //   this->previous_stress(el_type, ghost_type).begin(dim, dim);

  for (; strain_it != strain_end; ++strain_it, ++d_strain_it, ++sigma_it, ++sigma_th_cur_it, ++sigma_th_prev_it, ++inelas_strain_it, ++iso_hardening) {
    Matrix<Real> & grad_u = *strain_it;
    Matrix<Real> & grad_delta_u = *d_strain_it;
    Matrix<Real> & sigma_tensor = *sigma_it;
    Matrix<Real> & inelas_strain_tensor = *inelas_strain_it;
    //Real & iso_hard = *iso_hardening;


    //gradUToF<dim > (grad_u, F_tensor);
    /*switch (dim){
      case 3:
      Math::inv3(F_tensor.storage(), invF_tensor.storage());
      break;
      case 2:
      Math::inv2(F_tensor.storage(), invF_tensor.storage());
      break;

      }*/

    computeStressOnQuad(grad_u, grad_delta_u, sigma_tensor, inelas_strain_tensor,*iso_hardening, *sigma_th_cur_it, *sigma_th_prev_it);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialPlasticityinc<spatial_dimension>::computeTangentModuli(__attribute__((unused)) const ElementType & el_type,
								    Array<Real> & tangent_matrix,
								    __attribute__((unused)) GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Int dim=spatial_dimension;
  //     Array<Real>::iterator< Matrix<Real> > d_strain_it =
  //       this->delta_strain(el_type, ghost_type).begin(dim, dim);

  //  Array<Real>::iterator< Matrix<Real> > sigma_it =
  //  this->stress(el_type, ghost_type).begin(dim, dim);


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
Real MaterialPlasticityinc<spatial_dimension>::getPushWaveSpeed() const {
  return sqrt((lambda + 2 * mu) / this->rho);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialPlasticityinc<spatial_dimension>::getShearWaveSpeed() const {
  return sqrt(mu / this->rho);
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialPlasticityinc);

__END_AKANTU__
