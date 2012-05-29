/**
 * @file   material_elastic_caughey.cc
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
#include "material_elastic_caughey.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialElasticCaughey<spatial_dimension>::MaterialElasticCaughey(SolidMechanicsModel & model,
								  const ID & id)  :
  Material(model, id), MaterialElastic<spatial_dimension>(model, id),
  stress_viscosity("stress_viscosity", id),
  stress_elastic("stress_elastic", id) {
  AKANTU_DEBUG_IN();

  alpha = 0.;

  this->initInternalVector(this->stress_viscosity, spatial_dimension * spatial_dimension);
  this->initInternalVector(this->stress_elastic  , spatial_dimension * spatial_dimension);


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElasticCaughey<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialElastic<spatial_dimension>::initMaterial();

  this->resizeInternalVector(this->stress_viscosity);
  this->resizeInternalVector(this->stress_elastic  );

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElasticCaughey<spatial_dimension>::computeStress(ElementType el_type,
							      GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  Vector<UInt> & elem_filter = this->element_filter  (el_type, ghost_type);
  Vector<Real> & stress_visc = stress_viscosity(el_type, ghost_type);
  Vector<Real> & stress_el   = stress_elastic  (el_type, ghost_type);

  MaterialElastic<spatial_dimension>::computeStress(el_type, ghost_type);

  Vector<Real> & velocity = this->model->getVelocity();
  UInt nb_elem = this->element_filter(el_type, ghost_type).getSize();
  UInt nb_quad_points_per_elem =
    this->model->getFEM().getNbQuadraturePoints(el_type, ghost_type);
  UInt nb_quad_points = nb_quad_points_per_elem * nb_elem;

  Vector<Real> strain_rate(nb_quad_points,
			   spatial_dimension * spatial_dimension,
			   "strain_rate");

  this->model->getFEM().gradientOnQuadraturePoints(velocity, strain_rate,
   					     spatial_dimension,
   					     el_type, ghost_type, &elem_filter);

  Real F[3*3];
  Real sigma[3*3];
  Real * strain_rate_val = strain_rate.storage();
  Real * stress_visc_val = stress_visc.storage();
  Real * stress_el_val   = stress_el.storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;
  memset(F, 0, 3 * 3 * sizeof(Real));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      F[3*i + j] = strain_rate_val[spatial_dimension * i + j];

  MaterialElastic<spatial_dimension>::computeStress(F, sigma);

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j) {
      stress_visc_val[spatial_dimension * i + j]  = alpha * sigma[3 * i + j];
      stress_el_val  [spatial_dimension * i + j]  = stress_val[spatial_dimension*i + j];
      stress_val     [spatial_dimension * i + j] += stress_visc_val[spatial_dimension*i + j];
    }

  strain_rate_val += spatial_dimension * spatial_dimension;
  stress_visc_val += spatial_dimension * spatial_dimension;
  stress_el_val   += spatial_dimension * spatial_dimension;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElasticCaughey<spatial_dimension>::computePotentialEnergy(ElementType el_type,
								       GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if(ghost_type != _not_ghost) return;

  Vector<Real> & stress_el = stress_elastic(el_type, ghost_type);
  Real * stress_el_val = stress_el.storage();

  Real * epot = this->potential_energy(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;

  MaterialElastic<spatial_dimension>::computePotentialEnergy(strain_val, stress_el_val, epot);
  epot++;
  stress_el_val += spatial_dimension*spatial_dimension;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
bool MaterialElasticCaughey<spatial_dimension>::setParam(const std::string & key,
							 const std::string & value,
				      const ID & id) {
  std::stringstream sstr(value);
  if(key == "alpha") { sstr >> alpha; }
  else { return MaterialElastic<spatial_dimension>::setParam(key, value, id); }
  return true;
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElasticCaughey<spatial_dimension>::printself(std::ostream & stream,
							  int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "MaterialElasticCaughey [" << std::endl;
  MaterialElastic<spatial_dimension>::printself(stream, indent + 1);
  stream << space << " + artifical viscous ratio : " << alpha << std::endl;
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialElasticCaughey);


__END_AKANTU__
