/**
 * @file   material_mazars_non_local.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marion Chambart <marion.chambart@epfl.ch>
 * @date   Tue Jul 27 11:53:52 2010
 *
 * @brief  Specialization of the material class for the non-local mazars material
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
#include "material_mazars_non_local.hh"
#include "solid_mechanics_model.hh"


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialMazarsNonLocal::MaterialMazarsNonLocal(Model & model, const ID & id)  :
  Material(model, id), MaterialElastic(model, id),
  MaterialMazars(model, id), MaterialNonLocal(model, id),
  Ehat("Ehat", id) {
  AKANTU_DEBUG_IN();

  is_non_local = true;

  initInternalVector(this->Ehat, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialMazarsNonLocal::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialMazars::initMaterial();
  MaterialNonLocal::initMaterial();

  resizeInternalVector(this->Ehat);

  is_init = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialMazarsNonLocal::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real F[3*3];
  Real sigma[3*3];
  Real * dam = damage(el_type, ghost_type).storage();
  Real * Ehatt = Ehat(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;
  memset(F, 0, 3 * 3 * sizeof(Real));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      F[3*i + j] = strain_val[spatial_dimension * i + j];

  MaterialMazars::computeStress(F, sigma, *dam, *Ehatt);
  ++dam;
  ++Ehatt;

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      stress_val[spatial_dimension*i + j] = sigma[3 * i + j];

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialMazarsNonLocal::computeNonLocalStress(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  ByElementTypeReal Ehatnl("Ehat non local", id);
  initInternalVector(Ehatnl, 1);
  resizeInternalVector(Ehatnl);

  weigthedAvergageOnNeighbours(Ehat, Ehatnl, 1);

  UInt spatial_dimension = model->getSpatialDimension();

  Mesh::type_iterator it = model->getFEM().getMesh().firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = model->getFEM().getMesh().lastType(spatial_dimension, ghost_type);

  for(; it != last_type; ++it) {
    computeNonLocalStress(Ehatnl(*it, ghost_type), *it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialMazarsNonLocal::computeNonLocalStress(Vector<Real> & Ehatnl, ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam = damage(el_type, ghost_type).storage();
  Real * Ehatnlt = Ehatnl.storage();

  Real sigma[3*3];
  Real F[3*3] ;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;
  memset(F, 0, 3 * 3 * sizeof(Real));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j) {
      F[3*i + j] = strain_val[spatial_dimension * i + j];
      sigma[3 * i + j] = stress_val[spatial_dimension*i + j];
    }

  computeDamageAndStress(F, sigma, *dam, *Ehatnlt);

  ++dam;
  ++Ehatnlt;

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      stress_val[spatial_dimension*i + j] = sigma[3 * i + j];

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool MaterialMazarsNonLocal::setParam(const std::string & key, const std::string & value,
			       const ID & id) {


    return MaterialNonLocal::setParam(key, value, id) ||
      MaterialMazars::setParam(key, value, id);

}


/* -------------------------------------------------------------------------- */
void MaterialMazarsNonLocal::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_mazars_non_local> [" << std::endl;
  MaterialMazars::printself(stream, indent + 1);
  MaterialNonLocal::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
