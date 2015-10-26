/**
 * @file   fe_engine_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Jul 20 2010
 * @date last modification: Fri Jun 13 2014
 *
 * @brief  Implementation of the inline functions of the FEEngine Class
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline Real FEEngine::getElementInradius(const Matrix<Real> & coord, const ElementType & type) {
  AKANTU_DEBUG_IN();

  Real inradius = 0;

#define GET_INRADIUS(type)				\
  inradius = ElementClass<type>::getInradius(coord);	\

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_INRADIUS);
#undef GET_INRADIUS

  AKANTU_DEBUG_OUT();
  return inradius;
}

/* -------------------------------------------------------------------------- */
inline InterpolationType FEEngine::getInterpolationType(const ElementType & type) {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = _itp_not_defined;

#define GET_ITP(type)						\
  itp_type = ElementClassProperty<type>::interpolation_type;	\

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_ITP);
#undef GET_ITP

  AKANTU_DEBUG_OUT();
  return itp_type;
}

/* -------------------------------------------------------------------------- */
/// @todo rewrite this function in order to get the cohesive element
/// type directly from the facet
#if defined(AKANTU_COHESIVE_ELEMENT)
inline ElementType FEEngine::getCohesiveElementType(const ElementType & type_facet) {
  AKANTU_DEBUG_IN();

  ElementType type_cohesive = _not_defined;

  if (type_facet == _point_1) type_cohesive = _cohesive_1d_2;
  else if (type_facet == _segment_2) type_cohesive = _cohesive_2d_4;
  else if (type_facet == _segment_3) type_cohesive = _cohesive_2d_6;
  else if (type_facet == _triangle_3) type_cohesive = _cohesive_3d_6;
  else if (type_facet == _triangle_6) type_cohesive = _cohesive_3d_12;
  else if (type_facet == _quadrangle_4) type_cohesive = _cohesive_3d_8;
  else if (type_facet == _quadrangle_8) type_cohesive = _cohesive_3d_16;

  AKANTU_DEBUG_OUT();
  return type_cohesive;
}
#else
inline ElementType FEEngine::getCohesiveElementType(__attribute__((unused)) const ElementType & type_facet) {
  return _not_defined;
}
#endif

/* -------------------------------------------------------------------------- */
#if defined(AKANTU_IGFEM)
__END_AKANTU__
#include "igfem_helper.hh"
__BEGIN_AKANTU__

inline Vector<ElementType> FEEngine::getIGFEMElementTypes(const ElementType & type) {

#define GET_IGFEM_ELEMENT_TYPES(type) \
  return IGFEMHelper::getIGFEMElementTypes<type>();
  
 AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_IGFEM_ELEMENT_TYPES);

#undef GET_IGFEM_ELEMENT_TYPES
}
#endif

/* -------------------------------------------------------------------------- */
template<typename T>
void FEEngine::extractNodalToElementField(const Mesh & mesh,
					  const Array<T> & nodal_f,
					  Array<T> & elemental_f,
					  const ElementType & type,
					  const GhostType & ghost_type,
					  const Array<UInt> & filter_elements) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_degree_of_freedom = nodal_f.getNbComponent();
  UInt nb_element = mesh.getNbElement(type, ghost_type);
  UInt * conn_val = mesh.getConnectivity(type, ghost_type).storage();

  if(filter_elements != empty_filter) {
    nb_element      = filter_elements.getSize();
  }

  elemental_f.resize(nb_element);

  T * nodal_f_val = nodal_f.storage();
  T * f_val = elemental_f.storage();

  UInt * el_conn;
  for (UInt el = 0; el < nb_element; ++el) {
    if(filter_elements != empty_filter) el_conn = conn_val + filter_elements(el) * nb_nodes_per_element;
    else el_conn = conn_val + el * nb_nodes_per_element;

    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt node = *(el_conn + n);
      std::copy(nodal_f_val +      node  * nb_degree_of_freedom,
		nodal_f_val + (node + 1) * nb_degree_of_freedom,
		f_val);
      f_val += nb_degree_of_freedom;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<typename T>
void FEEngine::filterElementalData(const Mesh & mesh,
				   const Array<T> & elem_f,
				   Array<T> & filtered_f,
				   const ElementType & type,
				   const GhostType & ghost_type,
				   const Array<UInt> & filter_elements) {
  AKANTU_DEBUG_IN();

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  if(nb_element == 0) {
    filtered_f.resize(0);
    return;
  }

  UInt nb_degree_of_freedom = elem_f.getNbComponent();
  UInt nb_data_per_element  = elem_f.getSize() / nb_element;

  if(filter_elements != empty_filter) {
    nb_element      = filter_elements.getSize();
  }

  filtered_f.resize(nb_element * nb_data_per_element);

  T * elem_f_val = elem_f.storage();
  T * f_val = filtered_f.storage();

  UInt el_offset;
  for (UInt el = 0; el < nb_element; ++el) {
    if(filter_elements != empty_filter) el_offset = filter_elements(el);
    else el_offset = el;

    std::copy(elem_f_val +      el_offset  * nb_data_per_element * nb_degree_of_freedom,
	      elem_f_val + (el_offset + 1) * nb_data_per_element * nb_degree_of_freedom,
	      f_val);
    f_val += nb_degree_of_freedom * nb_data_per_element;
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
