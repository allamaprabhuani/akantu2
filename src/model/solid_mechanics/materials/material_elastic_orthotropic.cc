/**
 * @file   material_elastic_orthotropic.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Thu Apr 12 11:28:23 2012
 *
 * @brief  Orthotropic elastic material
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

/* -------------------------------------------------------------------------- */
#include "material_elastic_orthotropic.hh"
#include "solid_mechanics_model.hh"
#include <algorithm>

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialElasticOrthotropic::MaterialElasticOrthotropic(Model & model, const ID & id)  :
  Material(model, id) {
  AKANTU_DEBUG_IN();

  E1 = 0;
  E2 = 0;
  E3 = 0;
  nu12 = 0;
  nu13 = 0;
  nu23 = 0;
  G12 = 0;
  G13 = 0;
  G23 = 0;

  plane_stress = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialElasticOrthotropic::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();

  S = new Real[9];
  Real Delta;

  Real nu21 = nu12 * E2 / E1;
  Real nu31 = nu13 * E3 / E1;
  Real nu32 = nu23 * E3 / E2;

  if(plane_stress) {
    Delta = 1 - nu12 * nu21;

    S[0] = E1 / Delta;
    S[1] = E2 / Delta;
    S[2] = 0;
    S[3] = E1 * nu21 / Delta;
    S[4] = E2 * nu12 / Delta;
    S[5] = 0;
    S[6] = G12;
    S[7] = 0;
    S[8] = 0;
  }
  else {
    Delta = 1 - nu12 * nu21 - nu23 * nu32 - nu31 * nu13 - 2 * nu21 * nu13 * nu32;

    S[0] = E1 * (1 - nu23 * nu32) / Delta;
    S[1] = E2 * (1 - nu31 * nu13) / Delta;
    S[2] = E3 * (1 - nu12 * nu21) / Delta;

    S[3] = E2 * (nu12 + nu32 * nu13) / Delta;
    S[4] = E1 * (nu31 + nu21 * nu32) / Delta;
    S[5] = E3 * (nu23 + nu21 * nu13) / Delta;

    S[6] = G12;
    S[7] = G13;
    S[8] = G23;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialElasticOrthotropic::computeStress(ElementType el_type, GhostType ghost_type) {
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
void MaterialElasticOrthotropic::computeTangentStiffness(const ElementType & el_type,
					      Vector<Real> & tangent_matrix,
					      GhostType ghost_type) {

  switch(spatial_dimension) {
  case 1: { computeTangentStiffnessByDim<1>(el_type, tangent_matrix, ghost_type); break; }
  case 2: { computeTangentStiffnessByDim<2>(el_type, tangent_matrix, ghost_type); break; }
  case 3: { computeTangentStiffnessByDim<3>(el_type, tangent_matrix, ghost_type); break; }
  }
}

/* -------------------------------------------------------------------------- */
bool MaterialElasticOrthotropic::setParam(const std::string & key, const std::string & value,
			       const ID & id) {
  std::stringstream sstr(value);
  if(key == "E1") { sstr >> E1; }
  else if(key == "E2") { sstr >> E2; }
  else if(key == "E3") { sstr >> E3; }
  else if(key == "nu12") { sstr >> nu12; }
  else if(key == "nu13") { sstr >> nu13; }
  else if(key == "nu23") { sstr >> nu23; }
  else if(key == "G12") { sstr >> G12; }
  else if(key == "G13") { sstr >> G13; }
  else if(key == "G23") { sstr >> G23; }
  else if(key == "Plane_Stress") { sstr >> plane_stress; }
  else { return Material::setParam(key, value, id); }
  return true;
}

/* -------------------------------------------------------------------------- */
Real MaterialElasticOrthotropic::getPushWaveSpeed() {
  return sqrt( (*std::max_element(S, S+3)) /rho);
}

/* -------------------------------------------------------------------------- */
Real MaterialElasticOrthotropic::getShearWaveSpeed() {
  return sqrt( (*std::max_element(S+6, S+9)) /rho);
}

/* -------------------------------------------------------------------------- */
void MaterialElasticOrthotropic::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_elastic<_orthotropic> [" << std::endl;
  if(!plane_stress)
    stream << space << " + Plane strain" << std::endl;
  else
    stream << space << " + Plane stress" << std::endl;
  stream << space << " + density                 : " << rho << std::endl;
  stream << space << " + Young's modulus (x)     : " << E1 << std::endl;
  stream << space << " + Young's modulus (y)     : " << E2 << std::endl;
  stream << space << " + Young's modulus (z)     : " << E3 << std::endl;
  stream << space << " + Poisson's ratio (xy)    : " << nu12 << std::endl;
  stream << space << " + Poisson's ratio (xz)    : " << nu13 << std::endl;
  stream << space << " + Poisson's ratio (yz)    : " << nu23 << std::endl;
  stream << space << " + Shear modulus (xy)      : " << G12 << std::endl;
  stream << space << " + Shear modulus (xz)      : " << G13 << std::endl;
  stream << space << " + Shear modulus (yz)      : " << G23 << std::endl;

  Material::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

__END_AKANTU__
