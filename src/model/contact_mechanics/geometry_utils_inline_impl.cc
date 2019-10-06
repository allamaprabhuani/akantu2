#include "geometry_utils.hh"

#ifndef __AKANTU_GEOMETRY_UTILS_INLINE_IMPL_CC__
#define __AKANTU_GEOMETRY_UTILS_INLINE_IMPL_CC__


namespace akantu {

/* -------------------------------------------------------------------------- */
inline bool GeometryUtils::isBoundaryElement(const Mesh & mesh, const Element & element) {

  const auto & element_to_subelement =
    mesh.getElementToSubelement(element.type)(element.element);

  // for regular boundary elements
  if (element_to_subelement.size() == 1 and
      element_to_subelement[0].kind() == _ek_regular) {
    return true;
  }

  // for cohesive boundary elements
  UInt nb_subelements_regular = 0;
  for (auto subelem : element_to_subelement) {
    if (subelem == ElementNull)
      continue;
    
    if (subelem.kind() == _ek_regular)
      ++nb_subelements_regular;
  }

  auto nb_subelements = element_to_subelement.size();
  if (nb_subelements_regular < nb_subelements)
    return true;
     
  return false;
}

/* -------------------------------------------------------------------------- */
  inline bool GeometryUtils::isValidProjection(const Vector<Real> & projection) {
  
  UInt nb_xi_inside = 0;
  Real tolerance = 1e-3;

  for (auto xi : projection) {
    if (xi >= -1.0 - tolerance and xi <= 1.0 + tolerance)
      nb_xi_inside++;
  }

  if (nb_xi_inside == projection.size())
    return true;

  return false;
}

} //namespace akantu

#endif
