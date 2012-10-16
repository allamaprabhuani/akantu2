/**
 * @file   material_stiffness_proportional.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed May  4 15:25:52 2011
 *
 * @brief  Special. of the material class for the caughey viscoelastic material
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
#include "material_stiffness_proportional.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialStiffnessProportional<spatial_dimension>::MaterialStiffnessProportional(SolidMechanicsModel & model,
										const ID & id)  :
  Material(model, id), MaterialElastic<spatial_dimension>(model, id),
  stress_viscosity("stress_viscosity", id),
  stress_elastic("stress_elastic", id) {
  AKANTU_DEBUG_IN();

  this->registerParam("Alpha", alpha, 0., ParamAccessType(_pat_parsable | _pat_modifiable), "Artificial viscous ratio");

  this->initInternalVector(this->stress_viscosity, spatial_dimension * spatial_dimension);
  this->initInternalVector(this->stress_elastic  , spatial_dimension * spatial_dimension);


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialStiffnessProportional<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialElastic<spatial_dimension>::initMaterial();

  this->resizeInternalVector(this->stress_viscosity);
  this->resizeInternalVector(this->stress_elastic  );

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialStiffnessProportional<spatial_dimension>::computeStress(ElementType el_type,
							      GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  Vector<UInt> & elem_filter = this->element_filter  (el_type, ghost_type);
  Vector<Real> & stress_visc = stress_viscosity(el_type, ghost_type);
  Vector<Real> & stress_el   = stress_elastic  (el_type, ghost_type);

  MaterialElastic<spatial_dimension>::computeStress(el_type, ghost_type);

  Vector<Real> & velocity = this->model->getVelocity();

  Vector<Real> strain_rate(0, spatial_dimension * spatial_dimension,
			   "strain_rate");

  this->model->getFEM().gradientOnQuadraturePoints(velocity, strain_rate,
						   spatial_dimension,
						   el_type, ghost_type, &elem_filter);

  Vector<Real>::iterator<types::RMatrix> strain_rate_it =
    strain_rate.begin(spatial_dimension, spatial_dimension);
  Vector<Real>::iterator<types::RMatrix> stress_visc_it =
    stress_visc.begin(spatial_dimension, spatial_dimension);
  Vector<Real>::iterator<types::RMatrix> stress_el_it =
    stress_el.begin(spatial_dimension, spatial_dimension);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  types::RMatrix & grad_v     = *strain_rate_it;
  types::RMatrix & sigma_visc = *stress_visc_it;
  types::RMatrix & sigma_el   = *stress_el_it;

  MaterialElastic<spatial_dimension>::computeStressOnQuad(grad_v, sigma_visc);

  sigma_visc *= alpha;
  sigma_el.copy(sigma);
  sigma += sigma_visc;

  ++strain_rate_it;
  ++stress_visc_it;
  ++stress_el_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialStiffnessProportional<spatial_dimension>::computePotentialEnergy(ElementType el_type,
								       GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if(ghost_type != _not_ghost) return;

  Vector<Real> & stress_el = stress_elastic(el_type, ghost_type);
  Vector<Real>::iterator<types::RMatrix> stress_el_it =
    stress_el.begin(spatial_dimension, spatial_dimension);

  Real * epot = this->potential_energy(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  types::RMatrix & sigma_el    = *stress_el_it;
  Material::computePotentialEnergyOnQuad(grad_u,
					 sigma_el,
					 *epot);
  epot++;
  ++stress_el_it;
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */


INSTANSIATE_MATERIAL(MaterialStiffnessProportional);


__END_AKANTU__
