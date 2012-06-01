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
template<UInt spatial_dimension>
MaterialElastic<spatial_dimension>::MaterialElastic(SolidMechanicsModel & model, const ID & id)  :
  Material(model, id) {
  AKANTU_DEBUG_IN();

  E   = 0;
  nu  = 1./2.;
  plane_stress = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElastic<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();

  lambda   = nu * E / ((1 + nu) * (1 - 2*nu));
  mu       = E / (2 * (1 + nu));

  if(plane_stress) {
    //lambda = 2 * lambda * mu / (lambda + 2 * mu);
    lambda = nu * E / ((1 + nu)*(1 - nu));
  }

  kpa      = lambda + 2./3. * mu;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElastic<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
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
template<UInt spatial_dimension>
void MaterialElastic<spatial_dimension>::computeTangentStiffness(__attribute__((unused)) const ElementType & el_type,
								 Vector<Real> & tangent_matrix,
								 __attribute__((unused)) GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * tangent_val   = tangent_matrix.values;
  UInt offset_tangent  = tangent_matrix.getNbComponent();
  UInt nb_quads        = tangent_matrix.getSize();

  if (nb_quads == 0) return;

  memset(tangent_val, 0, offset_tangent * nb_quads * sizeof(Real));
  for (UInt q = 0; q < nb_quads; ++q, tangent_val += offset_tangent) {
    computeTangentStiffness(tangent_val);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
bool MaterialElastic<spatial_dimension>::setParam(const std::string & key, const std::string & value,
			       const ID & id) {
  std::stringstream sstr(value);
  if(key == "E") { sstr >> E; }
  else if(key == "nu") { sstr >> nu; }
  else if(key == "Plane_Stress") { sstr >> plane_stress; }
  else { return Material::setParam(key, value, id); }
  return true;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialElastic<spatial_dimension>::getParam(const ID & param) const {
  ID key = to_lower(param);
  if(key == "e") { return E; }
  else if(key == "nu") { return nu; }
  else if(key == "lambda") { return nu; }
  else if(key == "mu") { return mu; }
  else if(key == "kapa") { return kpa; }
  else return Material::getParam(param);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElastic<spatial_dimension>::setParam(const ID & param, Real value) {
  ID key = to_lower(param);
  if(key == "e") { E = value; }
  else if(key == "nu") { nu = value; }
  else if(key == "lambda") { nu = value; }
  else if(key == "mu") { mu = value; }
  else if(key == "kapa") { kpa = value; }
  else Material::setParam(param, value);
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialElastic<spatial_dimension>::getPushWaveSpeed() {
  return sqrt((lambda + 2*mu)/rho);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialElastic<spatial_dimension>::getShearWaveSpeed() {
  return sqrt(mu/rho);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialElastic<spatial_dimension>::printself(std::ostream & stream,
						   int indent) const {
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
  if(this->isInit()) {
    stream << space << " + First Lamé coefficient  : " << lambda << std::endl;
    stream << space << " + Second Lamé coefficient : " << mu << std::endl;
    stream << space << " + Bulk coefficient        : " << kpa << std::endl;
  }
  Material::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialElastic);

__END_AKANTU__
