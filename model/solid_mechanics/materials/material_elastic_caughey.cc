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
MaterialElasticCaughey::MaterialElasticCaughey(SolidMechanicsModel & model, const MaterialID & id)  :
  MaterialElastic(model, id) {
  AKANTU_DEBUG_IN();

  alpha = 0.;
  UInt spatial_dimension = model.getSpatialDimension();
  UInt stress_size = spatial_dimension * spatial_dimension;

  initInternalVector(this->ghost_stress_viscosity, stress_size, "stress_viscosity", _ghost);
  initInternalVector(this->stress_viscosity      , stress_size, "stress_viscosity");
  initInternalVector(this->ghost_stress_elastic  , stress_size, "stress_elastic", _ghost);
  initInternalVector(this->stress_elastic        , stress_size, "stress_elastic");


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialElasticCaughey::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialElastic::initMaterial();

  resizeInternalVector(this->ghost_stress_viscosity);
  resizeInternalVector(this->stress_viscosity      );  
  resizeInternalVector(this->ghost_stress_elastic  );
  resizeInternalVector(this->stress_elastic        );  

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialElasticCaughey::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  UInt spatial_dimension = model->getSpatialDimension();
  Vector<UInt> * elem_filter;
  Vector<Real> * stress_visc;
  Vector<Real> * stress_el  ;
  if(ghost_type == _not_ghost) {
    elem_filter = element_filter[el_type];
    stress_visc = stress_viscosity[el_type];
    stress_el   = stress_elastic  [el_type];
  } else {
    elem_filter = ghost_element_filter[el_type];
    stress_visc = ghost_stress_viscosity[el_type];
    stress_el   = ghost_stress_elastic  [el_type];      
  }


  MaterialElastic::computeStress(el_type, ghost_type);

  Vector<Real> & velocity = model->getVelocity();
  Vector<Real> strain_rate(0, spatial_dimension * spatial_dimension);
  model->getFEM().gradientOnQuadraturePoints(velocity, strain_rate,
   					     spatial_dimension,
   					     el_type, ghost_type, elem_filter);

  Real F[3*3];
  Real sigma[3*3];
  Real * strain_rate_val = strain_rate.values;
  Real * stress_visc_val = stress_visc->values;
  Real * stress_el_val   = stress_el->values;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;
  memset(F, 0, 3 * 3 * sizeof(Real));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      F[3*i + j] = strain_rate_val[spatial_dimension * i + j];

  MaterialElastic::computeStress(F, sigma);

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j) {
      stress_visc_val[spatial_dimension*i + j]  = alpha * sigma[3 * i + j];
      stress_el_val  [spatial_dimension*i + j]  = stress_val[spatial_dimension*i + j];
      stress_val     [spatial_dimension*i + j] += stress_visc_val[spatial_dimension*i + j];
    }
  strain_rate_val += spatial_dimension * spatial_dimension;
  stress_visc_val += spatial_dimension * spatial_dimension;
  stress_el_val   += spatial_dimension * spatial_dimension;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialElasticCaughey::computePotentialEnergy(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Vector<Real> * stress_el;
  if(ghost_type == _not_ghost) {
    stress_el = stress_elastic[el_type];
  } else {
    stress_el = ghost_stress_elastic[el_type];      
  }
  Real * stress_el_val = stress_el->values;

  if(ghost_type != _not_ghost) return;
  Real * epot = potential_energy[el_type]->values;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;

  MaterialElastic::computePotentialEnergy(strain_val, stress_el_val, epot);
  epot++;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialElasticCaughey::setParam(const std::string & key, const std::string & value,
			       const MaterialID & id) {
  std::stringstream sstr(value);
  if(key == "alpha") { sstr >> alpha; }
  else { MaterialElastic::setParam(key, value, id); }
}


/* -------------------------------------------------------------------------- */
void MaterialElasticCaughey::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "MaterialElasticCaughey [" << std::endl;
  MaterialElastic::printself(stream, indent + 1);
  stream << space << " + artifical viscous ratio : " << alpha << std::endl;
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

__END_AKANTU__
