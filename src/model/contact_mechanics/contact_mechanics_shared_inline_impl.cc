#ifndef __AKANTU_CONTACT_MECHANICS_SHARED_INLINE_IMPL_CC__
#define __AKANTU_CONTACT_MECHANICS_SHARED_INLINE_IMPL_CC__

/* -------------------------------------------------------------------------- */
#include "contact_mechanics_shared.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
inline Vector<Real> AbstractContactDetector::getNodePosition(UInt node) const {
  Vector<Real> position(spatial_dimension);
  for (UInt s : arange(spatial_dimension)) {
    position(s) = this->positions(node, s);
  }
  return position;
}

/* -------------------------------------------------------------------------- */
void AbstractContactDetector::coordinatesOfElement(
    const Element & el, Matrix<Real> & coords) const {
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);
  const Vector<UInt> connect = mesh.getConnectivity(el.type, _not_ghost)
                                   .begin(nb_nodes_per_element)[el.element];

  for (UInt n = 0; n < nb_nodes_per_element; ++n) {
    UInt node = connect[n];
    for (UInt s : arange(spatial_dimension)) {
      coords(s, n) = this->positions(node, s);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename vector_type>
inline auto AbstractContactDetector::computeElementSizes(vector_type && nodes) const {
  struct {
    Real min_size;
    Real max_size;
  } out_values;

  out_values = { std::numeric_limits<Real>::max(), std::numeric_limits<Real>::min() };

  for (auto node : nodes) {
    Array<Element> elements;
    this->mesh.getAssociatedElements(node, elements);

    for (auto element : elements) {
      UInt nb_nodes_per_element = mesh.getNbNodesPerElement(element.type);
      Matrix<Real> elem_coords(spatial_dimension, nb_nodes_per_element);
      this->coordinatesOfElement(element, elem_coords);

      Real elem_size = FEEngine::getElementInradius(elem_coords, element.type);
      out_values.max_size = std::max(out_values.max_size, elem_size);
      out_values.min_size = std::min(out_values.min_size, elem_size);
    }
  }

  return out_values;
}

} // namespace akantu

#endif // __AKANTU_CONTACT_MECHANICS_SHARED_INLINE_IMPL_CC__
