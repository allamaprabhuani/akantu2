/**
 * @file   material_elastic.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:53:52 2010
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
MaterialElastic::MaterialElastic(Model & model, const ID & id)  :
  Material(model, id) {
  AKANTU_DEBUG_IN();

  E   = 0;
  nu  = 1./2.;
  plane_stress = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialElastic::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();

  lambda   = nu * E / ((1 + nu) * (1 - 2*nu));
  mu       = E / (2 * (1 + nu));

  if(plane_stress)
    lambda = 2 * lambda * mu / (lambda + 2 * mu);

  kpa      = lambda + 2./3. * mu;


  is_init = true;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialElastic::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real F[3*3];
  Real sigma[3*3];

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;
  memset(F, 0, 3 * 3 * sizeof(Real));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      F[3*i + j] = strain_val[spatial_dimension * i + j];

  //  for (UInt i = 0; i < spatial_dimension; ++i) F[i*3 + i] -= 1;

  computeStress(F, sigma);

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      stress_val[spatial_dimension*i + j] = sigma[3 * i + j];

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialElastic::computeTangentStiffness(const ElementType & el_type,
					      Vector<Real> & tangent_matrix,
					      GhostType ghost_type) {

  switch(spatial_dimension) {
  case 1: { computeTangentStiffnessByDim<1>(el_type, tangent_matrix, ghost_type); break; }
  case 2: { computeTangentStiffnessByDim<2>(el_type, tangent_matrix, ghost_type); break; }
  case 3: { computeTangentStiffnessByDim<3>(el_type, tangent_matrix, ghost_type); break; }
  }
}

/* -------------------------------------------------------------------------- */
bool MaterialElastic::setParam(const std::string & key, const std::string & value,
			       const ID & id) {
  std::stringstream sstr(value);
  if(key == "E") { sstr >> E; }
  else if(key == "nu") { sstr >> nu; }
  else if(key == "Plane_Stress") { sstr >> plane_stress; }
  else { return Material::setParam(key, value, id); }
  return true;
}

/* -------------------------------------------------------------------------- */
Real MaterialElastic::getPushWaveSpeed() {
  return sqrt((lambda + 2*mu)/rho);
}

/* -------------------------------------------------------------------------- */
Real MaterialElastic::getShearWaveSpeed() {
  return sqrt(mu/rho);
}

/* -------------------------------------------------------------------------- */
void MaterialElastic::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_elastic> [" << std::endl;
  if(!plane_stress)
    stream << space << " + Plane strain" << std::endl;
  else
    stream << space << " + Plane stress" << std::endl;
  stream << space << " + density                 : " << rho << std::endl;
  stream << space << " + Young's modulus         : " << E << std::endl;
  stream << space << " + Poisson's ratio         : " << nu << std::endl;
  if(is_init) {
    stream << space << " + First Lamé coefficient  : " << lambda << std::endl;
    stream << space << " + Second Lamé coefficient : " << mu << std::endl;
    stream << space << " + Bulk coefficient        : " << kpa << std::endl;
  }
  Material::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

__END_AKANTU__
