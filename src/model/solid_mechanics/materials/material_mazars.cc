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
template<UInt spatial_dimension>
MaterialMazars<spatial_dimension>::MaterialMazars(SolidMechanicsModel & model,
						  const ID & id)  :
  Material(model, id),
  MaterialElastic<spatial_dimension>(model, id),
  MaterialDamage<spatial_dimension>(model, id),
  damage_in_compute_stress(true) {
  AKANTU_DEBUG_IN();
  K0   = 1e-4;
  At   = 0.8;
  Ac   = 1.4;
  Bc   = 1900;
  Bt   = 12000;
  beta = 1.06;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMazars<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialDamage<spatial_dimension>::initMaterial();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMazars<spatial_dimension>::computeStress(ElementType el_type,
						      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam = this->damage(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;

  Real Ehat;
  computeStressOnQuad(grad_u, sigma, *dam, Ehat);
  ++dam;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;


  if(!this->is_non_local) this->updateDissipatedEnergy(ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
bool MaterialMazars<spatial_dimension>::setParam(const std::string & key,
						 const std::string & value,
						 const ID & id) {
  std::stringstream sstr(value);
  if(key == "K0") { sstr >> K0; }
  else if(key == "At") { sstr >> At; }
  else if(key == "Bt") { sstr >> Bt; }
  else if(key == "Ac") { sstr >> Ac; }
  else if(key == "Bc") { sstr >> Bc; }
  else if(key == "beta") { sstr >> beta; }
  else { return MaterialDamage<spatial_dimension>::setParam(key, value, id); }
  return true;
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMazars<spatial_dimension>::printself(std::ostream & stream,
						  int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "MaterialMazars [" << std::endl;
  stream << space << " + K0    : " << K0 << std::endl;
  stream << space << " + At    : " << At << std::endl;
  stream << space << " + Bt    : " << Bt << std::endl;
  stream << space << " + Ac    : " << Ac << std::endl;
  stream << space << " + Bc    : " << Bc << std::endl;
  stream << space << " + beta  : " << beta << std::endl;
  MaterialDamage<spatial_dimension>::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialMazars);

__END_AKANTU__
