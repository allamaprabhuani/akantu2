/**
 * @file   material_damage_linear.cc
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
#include "material_damage_linear.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialDamageLinear::MaterialDamageLinear(Model & model, const ID & id)  :
  Material(model, id), MaterialElastic(model, id), MaterialDamage(model, id) {
  AKANTU_DEBUG_IN();

  Sigc = 1e5;
  Gc  = 2;

  is_non_local = false;

  initInternalVector(this->K, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialDamageLinear::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialDamage::initMaterial();

  resizeInternalVector(this->K);
  Epsmin = Sigc / E;
  Epsmax = 2 * Gc/ Sigc + Epsmin;

  const Mesh&mesh=model->getFEM().getMesh() ;
  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if (Mesh::getSpatialDimension(*it)!=mesh.getSpatialDimension() )
      continue ;
    UInt nb_element  = mesh.getNbElement(*it);
    UInt nb_quad = model->getFEM().getNbQuadraturePoints(*it);
    Vector <Real> & K_vec = K(*it);
    std::fill_n(K_vec.storage(),nb_element*nb_quad,Epsmin);
  }

  is_init = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialDamageLinear::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real F[3*3];
  Real sigma[3*3];
  Real * dam = damage(el_type, ghost_type).storage();
  Real * K = damage(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;
  memset(F, 0, 3 * 3 * sizeof(Real));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      F[3*i + j] = strain_val[spatial_dimension * i + j];

  computeStress(F, sigma,*dam,*K);
  ++dam;

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      stress_val[spatial_dimension*i + j] = sigma[3 * i + j];

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool MaterialDamageLinear::setParam(const std::string & key, const std::string & value,
			       const ID & id) {
  std::stringstream sstr(value);
  if(key == "Sigc") { sstr >> Sigc; }
  else if(key == "Gc") { sstr >> Gc; }
  else { return MaterialDamage::setParam(key, value, id); }
  return true;
}


/* -------------------------------------------------------------------------- */
void MaterialDamageLinear::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_damage_linear> [" << std::endl;
  stream << space << " + SigmaC : " << Sigc << std::endl;
  stream << space << " + GC     : " << Gc << std::endl;
  MaterialDamage::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

__END_AKANTU__
