/**
 * @file   material_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
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
#include <iostream>

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
inline void Material::gradUToF(const types::RMatrix & grad_u,
			       types::RMatrix & F) {
  UInt size_F = F.size();

  AKANTU_DEBUG_ASSERT(F.size() >= grad_u.size() && grad_u.size() == dim,
		      "The dimension of the tensor F should be greater or equal to the dimension of the tensor grad_u.");

  for (UInt i = 0; i < dim; ++i)
    for (UInt j = 0; j < dim; ++j)
      F(i, j) = grad_u(i, j);

  for (UInt i = 0; i < size_F; ++i) F(i, i) += 1;
}

/* -------------------------------------------------------------------------- */
inline void Material::rightCauchy(const types::RMatrix & F,
				  types::RMatrix & C) {
  C.mul<true, false>(F, F);
}

/* -------------------------------------------------------------------------- */
inline void Material::leftCauchy(const types::RMatrix & F,
				 types::RMatrix & B) {
  B.mul<false, true>(F, F);
}

/* -------------------------------------------------------------------------- */
inline void Material::computePotentialEnergyOnQuad(types::RMatrix & grad_u,
						   types::RMatrix & sigma,
						   Real & epot) {
  epot = 0.;
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      epot += sigma(i, j) * grad_u(i, j);

  epot *= .5;
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void Material::transferBMatrixToSymVoigtBMatrix(const types::RMatrix & B,
						       types::RMatrix & Bvoigt,
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
inline void Material::buildElementalFieldInterpolationCoodinates(__attribute__((unused)) const types::RMatrix & coordinates,
								 __attribute__((unused)) types::RMatrix & coordMatrix) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template<>
inline void Material::buildElementalFieldInterpolationCoodinates<_triangle_3>(const types::RMatrix & coordinates,
									      types::RMatrix & coordMatrix) {

  for (UInt i = 0; i < coordinates.rows(); ++i)
    coordMatrix(i, 0) = 1;
}

/* -------------------------------------------------------------------------- */
template<>
inline void Material::buildElementalFieldInterpolationCoodinates<_triangle_6>(const types::RMatrix & coordinates,
									      types::RMatrix & coordMatrix) {

  UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(_triangle_6);

  for (UInt i = 0; i < coordinates.rows(); ++i) {
    coordMatrix(i, 0) = 1;
    for (UInt j = 1; j < nb_quadrature_points; ++j)
      coordMatrix(i, j) = coordinates(i, j-1);
  }
}

/**
 * @todo Write a more efficient interpolation for quadrangles by
 * dropping unnecessary quadrature points
 *
 */

/* -------------------------------------------------------------------------- */
template<>
inline void Material::buildElementalFieldInterpolationCoodinates<_quadrangle_4>(const types::RMatrix & coordinates,
										types::RMatrix & coordMatrix) {

  for (UInt i = 0; i < coordinates.rows(); ++i) {
    Real x = coordinates(i, 0);
    Real y = coordinates(i, 1);

    coordMatrix(i, 0) = 1;
    coordMatrix(i, 1) = x;
    coordMatrix(i, 2) = y;
    coordMatrix(i, 3) = x * y;
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void Material::buildElementalFieldInterpolationCoodinates<_quadrangle_8>(const types::RMatrix & coordinates,
										types::RMatrix & coordMatrix) {

  for (UInt i = 0; i < coordinates.rows(); ++i) {

    UInt j = 0;
    Real x = coordinates(i, 0);
    Real y = coordinates(i, 1);

    for (UInt e = 0; e <= 2; ++e) {
      for (UInt n = 0; n <= 2; ++n) {
	coordMatrix(i, j) = std::pow(x, e) * std::pow(y, n);
	++j;
      }
    }

  }
}

/* -------------------------------------------------------------------------- */
template<ElementType type>
inline UInt Material::getSizeElementalFieldInterpolationCoodinates() {
  return model->getFEM().getNbQuadraturePoints(type);
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
inline UInt Material::getNbDataForElements(const Vector<Element> & elements,
					   SynchronizationTag tag) const {
  if(tag == _gst_smm_stress) {
    return spatial_dimension * spatial_dimension * sizeof(Real) * this->getNbQuadraturePoints(elements);
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
inline void Material::packElementData(CommunicationBuffer & buffer,
				      const Vector<Element> & elements,
				      SynchronizationTag tag) const {
  if(tag == _gst_smm_stress) {
    packElementDataHelper(stress, buffer, elements);
  }
}

/* -------------------------------------------------------------------------- */
inline void Material::unpackElementData(CommunicationBuffer & buffer,
					const Vector<Element> & elements,
					SynchronizationTag tag) {
  if(tag == _gst_smm_stress) {
    unpackElementDataHelper(stress, buffer, elements);
  }
}

/* -------------------------------------------------------------------------- */
inline UInt Material::getNbQuadraturePoints(const Vector<Element> & elements) const {
  UInt nb_quad = 0;
  Vector<Element>::const_iterator<Element> it  = elements.begin();
  Vector<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    nb_quad += this->model->getFEM().getNbQuadraturePoints(el.type, el.ghost_type);
  }
  return nb_quad;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void Material::packElementDataHelper(const ByElementTypeVector<T> & data_to_pack,
					    CommunicationBuffer & buffer,
					    const Vector<Element> & elements) const {
  packUnpackElementDataHelper<T, true>(const_cast<ByElementTypeVector<T> &>(data_to_pack),
				       buffer,
				       elements);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void Material::unpackElementDataHelper(ByElementTypeVector<T> & data_to_unpack,
					      CommunicationBuffer & buffer,
					      const Vector<Element> & elements) const {
  packUnpackElementDataHelper<T, false>(data_to_unpack, buffer, elements);
}

/* -------------------------------------------------------------------------- */
template<typename T, bool pack_helper>
inline void Material::packUnpackElementDataHelper(ByElementTypeVector<T> & data_to_pack,
						  CommunicationBuffer & buffer,
						  const Vector<Element> & element) const {
  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;
  UInt nb_quad_per_elem = 0;
  UInt nb_component = 0;

  Vector<T> * vect = NULL;
  Vector<UInt> * element_index_material = NULL;

  Vector<Element>::const_iterator<Element> it  = element.begin();
  Vector<Element>::const_iterator<Element> end = element.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    if(el.type != current_element_type || el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type   = el.ghost_type;

      vect = &data_to_pack(el.type, el.ghost_type);
      element_index_material = &(this->model->getElementIndexByMaterial(current_element_type, current_ghost_type));

      nb_quad_per_elem = this->model->getFEM().getNbQuadraturePoints(el.type, el.ghost_type);
      nb_component = vect->getNbComponent();
    }

    UInt el_id = (*element_index_material)(el.element, 0);
    types::Vector<T> data(vect->storage() + el_id * nb_component * nb_quad_per_elem,
			  nb_component * nb_quad_per_elem);
    if(pack_helper)
      buffer << data;
    else
      buffer >> data;
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
template<typename T>
void Material::removeQuadraturePointsFromVectors(ByElementTypeVector<T> & data,
						 const ByElementTypeUInt & new_numbering) {
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    ByElementTypeVector<UInt>::type_iterator it  = new_numbering.firstType(0, gt, _ek_not_defined);
    ByElementTypeVector<UInt>::type_iterator end = new_numbering.lastType(0, gt, _ek_not_defined);
    for (; it != end; ++it) {
      ElementType type = *it;
      if(data.exists(type, gt)){
	const Vector<UInt> & renumbering = new_numbering(type, gt);

	Vector<T> & vect = data(type, gt);

	UInt nb_quad_per_elem = this->model->getFEM().getNbQuadraturePoints(type, gt);
	UInt nb_component = vect.getNbComponent();

	Vector<T> tmp(renumbering.getSize()*nb_quad_per_elem, nb_component);
	UInt new_size = 0;
	for (UInt i = 0; i < vect.getSize(); ++i) {
	  UInt new_i = renumbering(i);
	  if(new_i != UInt(-1)) {
	    memcpy(tmp.storage() + new_i * nb_component * nb_quad_per_elem,
		   vect.storage() + i * nb_component * nb_quad_per_elem,
		   nb_component * nb_quad_per_elem * sizeof(T));
	    ++new_size;
	  }
	}
	tmp.resize(new_size * nb_quad_per_elem);
	vect.copy(tmp);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void Material::onElementsAdded(__attribute__((unused)) const Vector<Element> & element_list) {
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

/* -------------------------------------------------------------------------- */
inline void Material::onElementsRemoved(const Vector<Element> & element_list, const ByElementTypeUInt & new_numbering) {
  UInt my_num = model->getInternalIndexFromID(id);

  ByElementTypeUInt material_local_new_numbering("remove mat filter elem", id);

  Vector<Element>::const_iterator<Element> el_begin = element_list.begin();
  Vector<Element>::const_iterator<Element> el_end   = element_list.end();

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    ByElementTypeVector<UInt>::type_iterator it  = new_numbering.firstType(0, gt, _ek_not_defined);
    ByElementTypeVector<UInt>::type_iterator end = new_numbering.lastType(0, gt, _ek_not_defined);
    for (; it != end; ++it) {
      ElementType type = *it;
      if(element_filter.exists(type, gt)){
	Vector<UInt> & elem_filter = element_filter(type, gt);

	Vector<UInt> & element_index_material = model->getElementIndexByMaterial(type, gt);
	element_index_material.resize(model->getFEM().getMesh().getNbElement(type, gt)); // all materials will resize of the same size...

	if(!material_local_new_numbering.exists(type, gt))
	  material_local_new_numbering.alloc(elem_filter.getSize(), 1, type, gt);

	Vector<UInt> & mat_renumbering = material_local_new_numbering(type, gt);
	const Vector<UInt> & renumbering = new_numbering(type, gt);
	Vector<UInt> elem_filter_tmp;
	UInt ni = 0;
	Element el;
	el.type = type;
	el.ghost_type = gt;
	for (UInt i = 0; i < elem_filter.getSize(); ++i) {
	  el.element = elem_filter(i);
	  if(std::find(el_begin, el_end, el) == el_end) {
	    UInt new_el = renumbering(el.element);
	    AKANTU_DEBUG_ASSERT(new_el != UInt(-1), "A not removed element as been badly renumbered");
	    elem_filter_tmp.push_back(new_el);
	    mat_renumbering(i) = ni;
	    element_index_material(new_el, 0) = ni;
	    element_index_material(new_el, 1) = my_num;
	    ++ni;
	  } else {
	    mat_renumbering(i) = UInt(-1);
	  }
	}

	elem_filter.resize(elem_filter_tmp.getSize());
	elem_filter.copy(elem_filter);
      }
    }
  }

  for (std::map<ID, ByElementTypeReal *>::iterator it = internal_vectors_real.begin();
       it != internal_vectors_real.end();
       ++it) {
    this->removeQuadraturePointsFromVectors(*(it->second), material_local_new_numbering);
  }

  for (std::map<ID, ByElementTypeUInt *>::iterator it = internal_vectors_uint.begin();
       it != internal_vectors_uint.end();
       ++it) {
    this->removeQuadraturePointsFromVectors(*(it->second), material_local_new_numbering);
  }
}
