/**
 * @file   material_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:57:43 2010
 *
 * @brief  Implementation of the inline functions of the class material
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

__END_AKANTU__

#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
inline UInt Material::addElement(const ElementType & type,
				 UInt element,
				 const GhostType & ghost_type) {
  element_filter(type, ghost_type).push_back(element);
  return element_filter(type, ghost_type).getSize()-1;
}

/* -------------------------------------------------------------------------- */
inline UInt Material::getTangentStiffnessVoigtSize(UInt dim) const {
  return (dim * (dim - 1) / 2 + dim);
}
/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void Material::gradUToF(const types::Matrix & grad_u,
			       types::Matrix & F) {
  UInt size_F = F.size();

  AKANTU_DEBUG_ASSERT(F.size() >= grad_u.size() && grad_u.size() == dim,
		      "The dimension of the tensor F should be greater or equal to the dimension of the tensor grad_u.");

  for (UInt i = 0; i < dim; ++i)
    for (UInt j = 0; j < dim; ++j)
      F(i, j) = grad_u(i, j);

  for (UInt i = 0; i < size_F; ++i) F(i, i) += 1;
}

/* -------------------------------------------------------------------------- */
inline void Material::rightCauchy(const types::Matrix & F,
				  types::Matrix & C) {
  C.mul<true, false>(F, F);
}

/* -------------------------------------------------------------------------- */
inline void Material::leftCauchy(const types::Matrix & F,
				 types::Matrix & B) {
  B.mul<false, true>(F, F);
}

/* -------------------------------------------------------------------------- */
inline void Material::computePotentialEnergyOnQuad(types::Matrix & grad_u,
						   types::Matrix & sigma,
						   Real & epot) {
  epot = 0.;
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      epot += sigma(i, j) * grad_u(i, j);

  epot *= .5;
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void Material::transferBMatrixToSymVoigtBMatrix(const types::Matrix & B,
						       types::Matrix & Bvoigt,
						       UInt nb_nodes_per_element) const {
  Bvoigt.clear();

  for (UInt i = 0; i < dim; ++i)
    for (UInt n = 0; n < nb_nodes_per_element; ++n)
      Bvoigt(i, i + n*dim) = B(n, i);

  if(dim == 2) {
    ///in 2D, fill the @f$ [\frac{\partial N_i}{\partial x}, \frac{\partial N_i}{\partial y}]@f$ row
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      Bvoigt(2, 1 + n*2) = B(n, 0);
      Bvoigt(2, 0 + n*2) = B(n, 1);
    }
  }


  if(dim == 3) {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      Real dndx = B(n, 0);
      Real dndy = B(n, 1);
      Real dndz = B(n, 2);

      ///in 3D, fill the @f$ [0, \frac{\partial N_i}{\partial y}, \frac{N_i}{\partial z}]@f$ row
      Bvoigt(3, 1 + n*3) = dndz;
      Bvoigt(3, 2 + n*3) = dndy;

      ///in 3D, fill the @f$ [\frac{\partial N_i}{\partial x}, 0, \frac{N_i}{\partial z}]@f$ row
      Bvoigt(4, 0 + n*3) = dndz;
      Bvoigt(4, 2 + n*3) = dndx;

      ///in 3D, fill the @f$ [\frac{\partial N_i}{\partial x}, \frac{N_i}{\partial y}, 0]@f$ row
      Bvoigt(5, 0 + n*3) = dndy;
      Bvoigt(5, 1 + n*3) = dndx;
    }
  }
}

/* -------------------------------------------------------------------------- */
template<ElementType type>
inline void Material::buildElementalFieldInterpolationCoodinates(__attribute__((unused)) const types::Matrix & coordinates,
								 __attribute__((unused)) types::Matrix & coordMatrix) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template<>
inline void Material::buildElementalFieldInterpolationCoodinates<_triangle_3>(const types::Matrix & coordinates,
									      types::Matrix & coordMatrix) {

  for (UInt i = 0; i < coordinates.rows(); ++i)
    coordMatrix(i, 0) = 1;
}

/* -------------------------------------------------------------------------- */
template<>
inline void Material::buildElementalFieldInterpolationCoodinates<_triangle_6>(const types::Matrix & coordinates,
									      types::Matrix & coordMatrix) {

  UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(_triangle_6);

  for (UInt i = 0; i < coordinates.rows(); ++i) {
    coordMatrix(i, 0) = 1;
    for (UInt j = 1; j < nb_quadrature_points; ++j)
      coordMatrix(i, j) = coordinates(i, j-1);
  }
}

/* -------------------------------------------------------------------------- */
template<ElementType type>
inline UInt Material::getSizeElementalFieldInterpolationCoodinates() {
  return model->getFEM().getNbQuadraturePoints(type);
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
inline void Material::extractElementalFieldForInterplation(const Vector<Real> & field,
							   Vector<Real> & filtered_field) {
  filtered_field.copy(field);
}

/* -------------------------------------------------------------------------- */
template<typename T>
void Material::registerParam(std::string name, T & variable, T default_value,
			     ParamAccessType type,
			     std::string description) {
  AKANTU_DEBUG_IN();
  params.registerParam<T>(name, variable, default_value, type, description);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<typename T>
void Material::registerParam(std::string name, T & variable, ParamAccessType type,
			     std::string description) {
  AKANTU_DEBUG_IN();
  params.registerParam<T>(name, variable, type, description);
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
inline UInt Material::getNbDataToPack(const Element & element,
				      SynchronizationTag tag) const {
  UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(element.type);
  if(tag == _gst_smm_stress) return spatial_dimension * spatial_dimension * sizeof(Real) * nb_quadrature_points;
  return 0;
}

/* -------------------------------------------------------------------------- */
inline UInt Material::getNbDataToUnpack(const Element & element,
					SynchronizationTag tag) const {
  UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(element.type);
  if(tag == _gst_smm_stress) return spatial_dimension * spatial_dimension * sizeof(Real) * nb_quadrature_points;
  return 0;
}

/* -------------------------------------------------------------------------- */
inline void Material::packData(CommunicationBuffer & buffer,
			       const Element & element,
			       SynchronizationTag tag) const {
  if(tag == _gst_smm_stress) {
    UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(element.type);
    Vector<Real>::const_iterator<types::Matrix> stress_it = stress(element.type, _not_ghost).begin(spatial_dimension, spatial_dimension);
    stress_it += element.element * nb_quadrature_points;
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++stress_it)
      buffer << *stress_it;
  }
}

/* -------------------------------------------------------------------------- */
inline void Material::unpackData(CommunicationBuffer & buffer,
				 const Element & element,
				 SynchronizationTag tag) {
  if(tag == _gst_smm_stress) {
    UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(element.type);
    Vector<Real>::iterator<types::Matrix> stress_it = stress(element.type, _ghost).begin(spatial_dimension, spatial_dimension);
    stress_it += element.element * nb_quadrature_points;
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++stress_it)
      buffer >> *stress_it;
  }
}


/* -------------------------------------------------------------------------- */
template <typename T>
inline T Material::getParam(const ID & param) const {
  try {
    return params.get<T>(param);
  } catch (...) {
    AKANTU_EXCEPTION("No parameter " << param << " in the material " << id);
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void Material::setParam(const ID & param, T value) {
  try {
    params.set<T>(param, value);
  } catch(...) {
    AKANTU_EXCEPTION("No parameter " << param << " in the material " << id);
  }
  updateInternalParameters();
}

/* -------------------------------------------------------------------------- */
inline void Material::onElementsAdded(const Vector<Element> & element_list) {
  for (std::map<ID, ByElementTypeReal *>::iterator it = internal_vectors_real.begin();
       it != internal_vectors_real.end();
       ++it) {
    resizeInternalVector(*(it->second));
  }

  for (std::map<ID, ByElementTypeUInt *>::iterator it = internal_vectors_uint.begin();
       it != internal_vectors_uint.end();
       ++it) {
    resizeInternalVector(*(it->second));
  }
}
