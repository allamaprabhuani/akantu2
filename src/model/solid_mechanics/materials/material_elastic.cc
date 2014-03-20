/**
 * @file   material_elastic.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date   Tue Jul 27 18:15:37 2010
 *
 * @brief  Specialization of the material class for the elastic material
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
#include "material_elastic.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt dim>
MaterialElastic<dim>::MaterialElastic(SolidMechanicsModel & model, const ID & id)  :
  Material(model, id),
  MaterialThermal<dim>(model, id) {
  AKANTU_DEBUG_IN();

  this->registerParam("Plane_Stress",plane_stress, false, _pat_parsmod, "Is plane stress"        ); /// @todo Plane_Stress should not be possible to be modified after initMaterial (but before)
  this->registerParam("lambda"      ,lambda             , _pat_readable, "First Lamé coefficient" );
  this->registerParam("mu"          ,mu                 , _pat_readable, "Second Lamé coefficient");
  this->registerParam("kapa"        ,kpa                , _pat_readable, "Bulk coefficient"       );

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialElastic<dim>::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();
  MaterialThermal<dim>::initMaterial();

  if (dim == 1) this->nu = 0.;

  updateInternalParameters();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialElastic<dim>::updateInternalParameters() {
  MaterialThermal<dim>::updateInternalParameters();

  lambda   = this->nu * this->E / ((1 + this->nu) * (1 - 2*this->nu));
  mu       = this->E / (2 * (1 + this->nu));

  if(plane_stress) lambda = this->nu * this->E / ((1 + this->nu)*(1 - this->nu));

  kpa      = lambda + 2./3. * mu;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElastic<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MaterialThermal<spatial_dimension>::computeStress(el_type, ghost_type);
  Array<Real>::iterator<> sigma_th_it = this->sigma_th_cur(el_type, ghost_type).begin();

  if (!this->finite_deformation) {
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
    computeStressOnQuad(grad_u, sigma, *sigma_th_it);
    ++sigma_th_it;
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  } else {

    /// deformation gradient
    Matrix<Real> F(spatial_dimension, spatial_dimension);
    /// finite strains
    Matrix<Real> E(spatial_dimension, spatial_dimension);

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
    /// compute F
    this->template gradUToF<spatial_dimension> (grad_u, F);
    /// compute E
    this->gradUToGreenStrain(grad_u, E);

    /// compute second Piola-Kirchhoff stress tensor
    computeStressOnQuad(E, sigma, *sigma_th_it);
    /// convert it to its Cauchy stress tensor
    this->template computeCauchyStressOnQuad<spatial_dimension>(F, sigma, sigma);

    ++sigma_th_it;
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElastic<spatial_dimension>::computeTangentModuli(__attribute__((unused)) const ElementType & el_type,
							      Array<Real> & tangent_matrix,
							      __attribute__((unused)) GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  computeTangentModuliOnQuad(tangent);
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialElastic<spatial_dimension>::getPushWaveSpeed() const {
  return sqrt((lambda + 2*mu)/this->rho);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialElastic<spatial_dimension>::getShearWaveSpeed() const {
  return sqrt(mu/this->rho);
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialElastic);

__END_AKANTU__
