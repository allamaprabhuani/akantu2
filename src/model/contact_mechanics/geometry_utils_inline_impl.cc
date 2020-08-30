#include "geometry_utils.hh"

#ifndef __AKANTU_GEOMETRY_UTILS_INLINE_IMPL_CC__
#define __AKANTU_GEOMETRY_UTILS_INLINE_IMPL_CC__


namespace akantu {

/* -------------------------------------------------------------------------- */
inline bool GeometryUtils::isBoundaryElement(const Mesh & mesh,
					     const Element & subelement) {

  const auto & element_to_subelement =
    mesh.getElementToSubelement(subelement.type)(subelement.element);

  // for regular boundary elements when surfaceselector is set to
  // physical surfaces, the mesh contains only 1 element attached to a
  // boundary subelement 
  if (element_to_subelement.size() == 1 and
      element_to_subelement[0].kind() == _ek_regular) {
    return true;
  }

  // for cohesive interface elements when surfaceSelector is set
  // either cohesive surface selector or all surface selector, in this
  // case mesg passes is actually mesh_facet and for boundary or
  // cohesive  interface 2 elements are associated to a subelement
  // we want only one regular element attached to the subelement
  
  UInt nb_elements_regular = 0;
  UInt nb_elements_cohesive = 0;

  for (auto elem : element_to_subelement) {
    if (elem == ElementNull)
      continue;
    
    if (elem.kind() == _ek_regular)
      ++nb_elements_regular;

    if (elem.kind() == _ek_cohesive)
      ++nb_elements_cohesive;
  }

  auto nb_elements = element_to_subelement.size();
  if (nb_elements_regular  < nb_elements)
    return true;
     
  return false;
}

/* -------------------------------------------------------------------------- */
  inline bool GeometryUtils::isValidProjection(const Vector<Real> & projection, 
					       Real extension_tolerance) {
  
  UInt nb_xi_inside = 0;

  for (auto xi : projection) {
    if (xi >= -1.0 - extension_tolerance and xi <= 1.0 + extension_tolerance)
      nb_xi_inside++;
  }

  if (nb_xi_inside == projection.size())
    return true;

  return false;
}

} //namespace akantu

#endif
