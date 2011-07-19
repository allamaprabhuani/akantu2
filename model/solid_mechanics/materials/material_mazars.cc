/**
 * @file   material_mazars.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Chambart <marion.chambart@epfl.ch>
 * @date   Tue Jul 27 11:53:52 2010
 *
 * @brief  Specialization of the material class for the damage material
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
#include "material_mazars.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialMazars::MaterialMazars(Model & model, const ID & id)  :
  Material(model, id),
  damage("damage", id) {
  AKANTU_DEBUG_IN();

  E   = 0;
  nu  = 1./2.;
  K0  = 1e-4;
  At = 0.8;
  Ac = 1.4;
  Bc= 1900;
  Bt = 12000;
  beta = 1.06;

  initInternalVector(this->damage, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialMazars::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();

  resizeInternalVector(this->damage);

  is_init = true;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialMazars::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real F[3*3];
  Real sigma[3*3];
  Real * dam = damage(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;
  memset(F, 0, 3 * 3 * sizeof(Real));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      F[3*i + j] = strain_val[spatial_dimension * i + j];

  computeStress(F, sigma, *dam);
  ++dam;

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      stress_val[spatial_dimension*i + j] = sigma[3 * i + j];

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialMazars::setParam(const std::string & key, const std::string & value,
			      const ID & id) {
  std::stringstream sstr(value);
  if(key == "E") { sstr >> E; }
  else if(key == "nu") { sstr >> nu; }
  else if(key == "K0") { sstr >> K0; }
  else if(key == "At") { sstr >> At; }
  else if(key == "Bt") { sstr >> Bt; }
  else if(key == "Ac") { sstr >> Ac; }
  else if(key == "Bc") { sstr >> Bc; }
  else if(key == "beta") { sstr >> beta; }
  else { Material::setParam(key, value, id); }
}


/* -------------------------------------------------------------------------- */
void MaterialMazars::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_mazars> [" << std::endl;
  stream << space << " + id                      : " << id << std::endl;
  stream << space << " + name                    : " << name << std::endl;
  stream << space << " + density                 : " << rho << std::endl;
  stream << space << " + Young's modulus         : " << E << std::endl;
  stream << space << " + Poisson's ratio         : " << nu << std::endl;
  stream << space << " + K0                      : " << K0 << std::endl;
  stream << space << " + At                      : " << At << std::endl;
  stream << space << " + Bt                      : " << Bt << std::endl;
  stream << space << " + Ac                      : " << Ac << std::endl;
  stream << space << " + Bc                      : " << Bc << std::endl;
  stream << space << " + beta                    : " << beta << std::endl;
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

__END_AKANTU__
