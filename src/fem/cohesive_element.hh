/**
 * @file   cohesive_element.hh
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Feb  2 14:07:49 2012
 *
 * @brief  Cohesive element class
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
#include "element_class.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COHESIVE_ELEMENT_HH__
#define __AKANTU_COHESIVE_ELEMENT_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<ElementType t>
struct CohesiveElementSubElementType {
  enum { value = _not_defined };
};

template<>
struct CohesiveElementSubElementType<_cohesive_2d_4> {
  enum { value = _segment_2 };
};

template<>
struct CohesiveElementSubElementType<_cohesive_2d_6> {
  enum { value = _segment_3 };
};


/* -------------------------------------------------------------------------- */
template<ElementType ct>
class CohesiveElement : public ElementClass<ElementType(CohesiveElementSubElementType<ct>::value)> {
  typedef ElementClass<ElementType(CohesiveElementSubElementType<ct>::value)> Parent;

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {};

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  static UInt getNbNodesPerElement() {
    return 2 * Parent::getNbNodesPerElement();
  }

  static UInt getNbQuadraturePoints() {
    return Parent::getNbQuadraturePoints();
  }

  static UInt getSpatialDimension() {
    return Parent::getSpatialDimension() + 1;
  }

  static UInt getNbShapeFunctions() {
    return Parent::getNbShapeFunctions();
  }

  static inline Real * getQuadraturePoints() {
    return Parent::getQuadraturePoints();
  }

  static inline UInt getShapeSize() {
    return Parent::getShapeSize();
  }

  static inline UInt getShapeDerivativesSize(){
    return Parent::getShapeDerivativesSize();
  }

  /// compute the in-radius
  static inline Real getInradius(const Real * coord) {
    return Parent::getInradius(coord);
  }

  static inline Real * getGaussIntegrationWeights() {
    return Parent::getGaussIntegrationWeights();
  }

  static inline ElementType getFacetElementType() {
    return ElementType(CohesiveElementSubElementType<ct>::value);
  }

  static inline UInt getNbFacetsPerElement() { return 2; }

  static AKANTU_GET_MACRO_NOT_CONST(FacetLocalConnectivityPerElement, facet_connectivity, UInt**);


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// local connectivity of facets
  static UInt * facet_connectivity[];

  /// vectorial connectivity of facets
  static UInt vec_facet_connectivity[];

};


__END_AKANTU__

#endif /* __AKANTU_COHESIVE_ELEMENT_HH__ */
