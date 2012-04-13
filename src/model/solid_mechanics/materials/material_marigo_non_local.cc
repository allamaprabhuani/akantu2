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
  Y("Y", id), update_weigths(0), compute_stress_calls(0) {
  AKANTU_DEBUG_IN();

  is_non_local = true;

  initInternalVector(this->Y, 1);

  AKANTU_DEBUG_OUT();
}

// /* -------------------------------------------------------------------------- */
// template<>
// void MaterialMarigoNonLocal::initWeightFuncion<BaseWeightFunction>() {
//   BaseWeightFunction * weight_function = new BaseWeightFunction(radius);
//   MaterialNonLocalParent::initMaterial(*weight_function);
// }

// /* -------------------------------------------------------------------------- */
// template<>
// void MaterialMarigoNonLocal::initWeightFuncion<DamagedWeightFunction>() {
//   DamagedWeightFunction * weight_function = new DamagedWeightFunction(radius, damage);
//   MaterialNonLocalParent::initMaterial(*weight_function);
// }

/* -------------------------------------------------------------------------- */
template<>
void MaterialMarigoNonLocal::initWeightFuncion<StressBasedWeightFunction>() {
  StressBasedWeightFunction * weight_function = new StressBasedWeightFunction(radius,
									      10e9,
									      spatial_dimension,
									      *this);

  const Mesh & mesh = model->getFEM().getMesh();
  ByElementTypeReal quadrature_points_coordinates("quadrature_points_coordinates", id);
  mesh.initByElementTypeVector(quadrature_points_coordinates, spatial_dimension, 0);

  computeQuadraturePointsCoordinates(mesh.getNodes(), quadrature_points_coordinates);
  weight_function->setQuadraturePointsCoordinates(quadrature_points_coordinates);
  weight_function->updatePrincipalStress(_not_ghost);
  weight_function->updatePrincipalStress(_ghost);

  MaterialNonLocalParent::initMaterial(*weight_function);
}

/* -------------------------------------------------------------------------- */
void MaterialMarigoNonLocal::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialMarigo::initMaterial();

  initWeightFuncion<MarigoNonLocalWeightFunction>();

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
  Real * dpe = dissipated_energy(el_type, ghost_type).storage();

  Real delta_t = model->getTimeStep();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;
  memset(F, 0, 3 * 3 * sizeof(Real));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      F[3*i + j] = strain_val[spatial_dimension * i + j];

  MaterialMarigo::computeStress(F, sigma, *dam, *Yt, *Ydq, delta_t, *dpe);
  ++dam;
  ++Yt;
  ++Ydq;
  ++dpe;

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      stress_val[spatial_dimension*i + j] = sigma[3 * i + j];

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialMarigoNonLocal::computeNonLocalStress(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const Mesh & mesh = model->getFEM().getMesh();

  if(update_weigths && (compute_stress_calls % update_weigths == 0)) {
    ByElementTypeReal quadrature_points_coordinates("quadrature_points_coordinates", id);
    mesh.initByElementTypeVector(quadrature_points_coordinates, spatial_dimension, 0);
    Vector<Real> coordinates(mesh.getNodes(), true);
    coordinates += model->getDisplacement();

    if(StressBasedWeightFunction * wf =
       dynamic_cast<StressBasedWeightFunction *>(weight_func)){
      computeQuadraturePointsCoordinates(coordinates, quadrature_points_coordinates);
      wf->setQuadraturePointsCoordinates(quadrature_points_coordinates);
      wf->updatePrincipalStress(ghost_type);
    }
    computeWeights(quadrature_points_coordinates);
  }

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
  Vector<Real>::iterator<types::Matrix> strain_it =
    strain(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

  Real * dam = damage(el_type, ghost_type).storage();
  Real * Ynlt = Ynl.storage();
  Real * Ydq = Yd_rand(el_type, ghost_type).storage();
  Real * dpe = dissipated_energy(el_type, ghost_type).storage();

  Real delta_t = model->getTimeStep();

  Real sigma[3*3];
  Real F[3*3];

  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      for (UInt i = 0; i < spatial_dimension; ++i)
        for (UInt j = 0; j < spatial_dimension; ++j) {
	  F[3*i + j]       = (*strain_it)(i, j);
	  sigma[3 * i + j] = (*stress_it)(i, j);
	}

      computeDamageAndStress(F, sigma, *dam, *Ynlt, *Ydq, delta_t, *dpe);

      for (UInt i = 0; i < spatial_dimension; ++i)
        for (UInt j = 0; j < spatial_dimension; ++j)
          (*stress_it)(i, j) = sigma[3 * i + j];

      ++stress_it;
      ++dam;
      ++Ynlt;
      ++Ydq;
      ++dpe;
    }
  }

  updateDissipatedEnergy(ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool MaterialMarigoNonLocal::setParam(const std::string & key, const std::string & value,
                               const ID & id) {
  std::stringstream sstr(value);
  if(key == "UpdateWeights") { sstr >> update_weigths; }
  else {
    return MaterialNonLocalParent::setParam(key, value, id) ||
      MaterialMarigo::setParam(key, value, id);
  }
  return true;
}


/* -------------------------------------------------------------------------- */
void MaterialMarigoNonLocal::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_marigo_non_local> [" << std::endl;
  stream << space << " + UpdateWeights  : " << update_weigths << std::endl;
  MaterialMarigo::printself(stream, indent + 1);
  MaterialNonLocalParent::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

__END_AKANTU__
