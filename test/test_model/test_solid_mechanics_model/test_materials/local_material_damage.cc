/**
 * @file   local_material_damage.cc
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
#include "local_material_damage.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
LocalMaterialDamage::LocalMaterialDamage(SolidMechanicsModel & model,
					 const ID & id)  :
  Material(model, id),
  damage("damage", id) {
  AKANTU_DEBUG_IN();

  E   = 0;
  nu  = 1./2.;
  Yd  = 50;
  Sd  = 5000;

  initInternalVector(this->damage, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void LocalMaterialDamage::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();

  resizeInternalVector(this->damage);

  lambda = nu * E / ((1 + nu) * (1 - 2*nu));
  mu     = E / (2 * (1 + nu));
  kpa    = lambda + 2./3. * mu;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void LocalMaterialDamage::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  resizeInternalVector(this->damage);

  Real * dam = damage(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;

  computeStressOnQuad(grad_u, sigma, *dam);
  ++dam;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void LocalMaterialDamage::computePotentialEnergy(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if(ghost_type != _not_ghost) return;
  Real * epot = potential_energy(el_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;

  computePotentialEnergyOnQuad(grad_u, sigma, *epot);
  epot++;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
bool LocalMaterialDamage::setParam(const std::string & key, const std::string & value,
			       const ID & id) {
  std::stringstream sstr(value);

  if(key == "E") { sstr >> E; }
  else if(key == "nu") { sstr >> nu; }
  else if(key == "Yd") { sstr >> Yd; }
  else if(key == "Sd") { sstr >> Sd; }
  else { return Material::setParam(key, value, id); }
  return true;
}


/* -------------------------------------------------------------------------- */
void LocalMaterialDamage::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_damage> [" << std::endl;
  stream << space << " + id                      : " << id << std::endl;
  stream << space << " + name                    : " << name << std::endl;
  stream << space << " + density                 : " << rho << std::endl;
  stream << space << " + Young's modulus         : " << E << std::endl;
  stream << space << " + Poisson's ratio         : " << nu << std::endl;
  stream << space << " + Yd                      : " << Yd << std::endl;
  stream << space << " + Sd                      : " << Sd << std::endl;
  if(this->isInit()) {
    stream << space << " + First Lamé coefficient  : " << lambda << std::endl;
    stream << space << " + Second Lamé coefficient : " << mu << std::endl;
    stream << space << " + Bulk coefficient        : " << kpa << std::endl;
  }
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

__END_AKANTU__
