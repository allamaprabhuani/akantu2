/**
 * @file   material_damage.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
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
#include "material_damage.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialDamage::MaterialDamage(SolidMechanicsModel & model, const MaterialID & id)  :
  Material(model, id) {
  AKANTU_DEBUG_IN();

  rho = 0;
  E   = 0;
  nu  = 1./2.;
  Yd  = 50;
  Sd  = 5000;

  for(UInt t = _not_defined; t < _max_element_type; ++t)
    this->damage[t] = NULL;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialDamage::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();

  const Mesh::ConnectivityTypeList & type_list =
    model->getFEM().getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;
    std::stringstream sstr_damage; sstr_damage << id << ":damage:" << *it;
    damage[*it] = &(alloc<Real>(sstr_damage.str(), 0,
				1, REAL_INIT_VALUE));
  }

  lambda = nu * E / ((1 + nu) * (1 - 2*nu));
  mu     = E / (2 * (1 + nu));
  kpa    = lambda + 2./3. * mu;

  is_init = true;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialDamage::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real F[3*3];
  Real sigma[3*3];
  damage[el_type]->resize(FEM::getNbQuadraturePoints(el_type)*element_filter[el_type]->getSize());
  Real * dam = damage[el_type]->values;

  MATERIAL_QUADRATURE_POINT_LOOP_BEGIN;
  memset(F, 0, 3 * 3 * sizeof(Real));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      F[3*i + j] = strain_val[spatial_dimension * i + j];

  for (UInt i = 0; i < spatial_dimension; ++i) F[i*3 + i] -= 1;

  computeStress(F, sigma,*dam);
  ++dam;

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      stress_val[spatial_dimension*i + j] = sigma[3 * i + j];

  MATERIAL_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialDamage::computePotentialEnergy(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if(ghost_type != _not_ghost) return;
  Real * epot = potential_energy[el_type]->values;

  MATERIAL_QUADRATURE_POINT_LOOP_BEGIN;

  computePotentialEnergy(strain_val, stress_val, epot);
  epot++;

  MATERIAL_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void MaterialDamage::setParam(const std::string & key, const std::string & value,
			       const MaterialID & id) {
  std::stringstream sstr(value);
  if(key == "rho") { sstr >> rho; }
  else if(key == "E") { sstr >> E; }
  else if(key == "nu") { sstr >> nu; }
  else if(key == "Yd") { sstr >> Yd; }
  else if(key == "Sd") { sstr >> Sd; }
  else { Material::setParam(key, value, id); }
}


/* -------------------------------------------------------------------------- */
void MaterialDamage::printself(std::ostream & stream, int indent) const {
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
  if(is_init) {
    stream << space << " + First Lamé coefficient  : " << lambda << std::endl;
    stream << space << " + Second Lamé coefficient : " << mu << std::endl;
    stream << space << " + Bulk coefficient        : " << kpa << std::endl;
  }
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

__END_AKANTU__
