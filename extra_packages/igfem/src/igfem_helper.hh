/**
 * @file   dumper_igfem_element_iterator.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Helper class to return sub element information
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

#ifndef __AKANTU_IGFEM_HELPER_HH__
#define __AKANTU_IGFEM_HELPER_HH__
/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

class FEEngine;

template <ElementType type> struct ElementTypeIGFEMData {
  static ElementType igfem_element_types[];
  static UInt nb_igfem_types;
};

struct IGFEMHelper {
  template <ElementType type>
  static VectorProxy<ElementType> getIGFEMElementTypes() {
    return VectorProxy<ElementType>(
        ElementTypeIGFEMData<type>::igfem_element_types,
        ElementTypeIGFEMData<type>::nb_igfem_types);
  }

  /// get the number of nodes for a given sub-element
  static UInt getNbNodesPerSubElement(const ElementType & type,
                                      const UInt sub_element) {
    UInt nb_nodes_per_sub_element = 0;
#define GET_NB_NODES_PER_SUB_ELEMENT(type)                                     \
  switch (sub_element) {                                                       \
  case 0:                                                                      \
    nb_nodes_per_sub_element =                                                 \
        ElementClass<ElementClassProperty<type>::sub_element_type_1>::         \
            getNbNodesPerInterpolationElement();                               \
    break;                                                                     \
  case 1:                                                                      \
    nb_nodes_per_sub_element =                                                 \
        ElementClass<ElementClassProperty<type>::sub_element_type_2>::         \
            getNbNodesPerInterpolationElement();                               \
    break;                                                                     \
  }
    AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(GET_NB_NODES_PER_SUB_ELEMENT);
#undef GET_NB_NODES_PER_SUB_ELEMENT
    return nb_nodes_per_sub_element;
  }

  /// get the connectivity for a given sub-element
  static UInt * getSubElementConnectivity(const ElementType & type,
                                          const UInt sub_element) {
    UInt * sub_element_connectivity = NULL;
#define GET_SUB_ELEMENT_CONNECTIVITY(type)                                     \
  sub_element_connectivity =                                                   \
      ElementClass<type>::getSubElementConnectivity(sub_element);
    AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(GET_SUB_ELEMENT_CONNECTIVITY);
#undef GET_SUB_ELEMENT_CONNECTIVITY
    return sub_element_connectivity;
  }

  /// get the sub-element type
  static ElementType getSubElementType(const ElementType & type,
                                       const UInt sub_element) {
    ElementType sub_type = _not_defined;
#define GET_SUB_ELEMENT_TYPE(type)                                             \
  switch (sub_element) {                                                       \
  case 0:                                                                      \
    sub_type = ElementClassProperty<type>::sub_element_type_1;                 \
    break;                                                                     \
  case 1:                                                                      \
    sub_type = ElementClassProperty<type>::sub_element_type_2;                 \
    break;                                                                     \
  }
    AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(GET_SUB_ELEMENT_TYPE);
#undef GET_SUB_ELEMENT_TYPE
    return sub_type;
  }

  /// get the nb of quads for one sub element type
  static UInt getNbQuadraturePoints(const ElementType & type,
                                    const UInt sub_element) {
    UInt nb_quad_points = 0;
#define GET_NB_QUADS(type)                                                     \
  switch (sub_element) {                                                       \
  case 0:                                                                      \
    nb_quad_points = GaussIntegrationElement<ElementClassProperty<             \
        type>::sub_element_type_1>::getNbQuadraturePoints();                   \
    break;                                                                     \
  case 1:                                                                      \
    nb_quad_points = GaussIntegrationElement<ElementClassProperty<             \
        type>::sub_element_type_2>::getNbQuadraturePoints();                   \
    break;                                                                     \
  }
    AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(GET_NB_QUADS);
#undef GET_NB_QUADS
    return nb_quad_points;
  }

  /// get the nb of parent nodes of a given igfem element type
  static UInt getNbParentNodes(const ElementType & type) {
    UInt nb_parent_nodes = 0;
#define GET_NB_PARENT_NODES(type)                                              \
  nb_parent_nodes =                                                            \
      ElementClass<ElementClassProperty<type>::parent_element_type>::          \
          getNbNodesPerInterpolationElement();

    AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(GET_NB_PARENT_NODES);
#undef GET_NB_PARENT_NODES
    return nb_parent_nodes;
  }

  /// get the nb of parent nodes of a given igfem element type
  static UInt getNbEnrichedNodes(const ElementType & type) {
    UInt nb_enriched_nodes = 0;
#define GET_NB_ENRICHED_NODES(type)                                            \
  nb_enriched_nodes = ElementClass<type>::getNbEnrichments();

    AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(GET_NB_ENRICHED_NODES);
#undef GET_NB_ENRICHED_NODES
    return nb_enriched_nodes;
  }

  /// get the nb of quads for one sub element type
  static UInt getElementOrientation(const ElementType & type,
                                    const Vector<bool> & is_inside) {
    UInt orientation = 0;

#define GET_ORIENTATION(type)                                                  \
  orientation = ElementClass<type>::getOrientation(is_inside);

    AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(GET_ORIENTATION);
#undef GET_ORIENTATION
    return orientation;
  }
};

__END_AKANTU__
#endif /* __AKANTU_IGFEM_HELPER_HH__ */
