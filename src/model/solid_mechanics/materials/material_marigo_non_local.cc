/**
 * @file   material_marigo.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Chambart <marion.chambart@epfl.ch>
 * @date   Tue Jul 27 11:53:52 2010
 *
 * @brief  Specialization of the material class for the marigo material
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
#include "material_marigo_non_local.hh"
#include "solid_mechanics_model.hh"


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialMarigoNonLocal::MaterialMarigoNonLocal(Model & model, const ID & id)  :
  Material(model, id), MaterialElastic(model, id),
  MaterialMarigo(model, id), MaterialNonLocalParent(model, id),
  Y("Y", id) {
  AKANTU_DEBUG_IN();

  is_non_local = true;

  initInternalVector(this->Y, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialMarigoNonLocal::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialMarigo::initMaterial();


  if(StressBasedWeightFunction * wf =
     dynamic_cast<StressBasedWeightFunction *>(weight_func)){
    wf->updatePrincipalStress(_not_ghost);
    wf->updatePrincipalStress(_ghost);
  }

  MaterialNonLocalParent::initMaterial();

  resizeInternalVector(this->Y);

  is_init = true;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void MaterialMarigoNonLocal::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real F[3*3];
  Real sigma[3*3];
  Real * dam = damage(el_type, ghost_type).storage();
  Real * Yt = Y(el_type, ghost_type).storage();
  Real * Ydq = Yd_rand(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;
  memset(F, 0, 3 * 3 * sizeof(Real));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      F[3*i + j] = strain_val[spatial_dimension * i + j];

  MaterialMarigo::computeStress(F, sigma, *dam, *Yt, *Ydq);
  ++dam;
  ++Yt;
  ++Ydq;

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      stress_val[spatial_dimension*i + j] = sigma[3 * i + j];

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialMarigoNonLocal::computeNonLocalStress(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  ByElementTypeReal Ynl("Y non local", id);
  initInternalVector(Ynl, 1);
  resizeInternalVector(Ynl);

  weightedAvergageOnNeighbours(Y, Ynl, 1);

  UInt spatial_dimension = model->getSpatialDimension();

  Mesh::type_iterator it = model->getFEM().getMesh().firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = model->getFEM().getMesh().lastType(spatial_dimension, ghost_type);

  for(; it != last_type; ++it) {
    computeNonLocalStress(Ynl(*it, ghost_type), *it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialMarigoNonLocal::computeNonLocalStress(Vector<Real> & Ynl,
                                                   ElementType el_type,
                                                   GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(el_type, ghost_type);
  UInt nb_element = element_filter(el_type, ghost_type).getSize();
  if (nb_element == 0) return;

  Vector<Real>::iterator<types::Matrix> stress_it =
    stress(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

  Real * dam = damage(el_type, ghost_type).storage();
  Real * Ynlt = Ynl.storage();
  Real * Ydq = Yd_rand(el_type, ghost_type).storage();

  Real sigma[3*3];

  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      for (UInt i = 0; i < spatial_dimension; ++i)
        for (UInt j = 0; j < spatial_dimension; ++j) {
	  sigma[3 * i + j] = (*stress_it)(i, j);
	}

      computeDamageAndStress(sigma, *dam, *Ynlt, *Ydq);

      for (UInt i = 0; i < spatial_dimension; ++i)
        for (UInt j = 0; j < spatial_dimension; ++j)
          (*stress_it)(i, j) = sigma[3 * i + j];

      ++stress_it;
      ++dam;
      ++Ynlt;
      ++Ydq;
    }
  }

  updateDissipatedEnergy(ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool MaterialMarigoNonLocal::setParam(const std::string & key, const std::string & value,
                               const ID & id) {
  return MaterialNonLocalParent::setParam(key, value, id) ||
    MaterialMarigo::setParam(key, value, id);
}


/* -------------------------------------------------------------------------- */
void MaterialMarigoNonLocal::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_marigo_non_local> [" << std::endl;
  MaterialMarigo::printself(stream, indent + 1);
  MaterialNonLocalParent::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

__END_AKANTU__
