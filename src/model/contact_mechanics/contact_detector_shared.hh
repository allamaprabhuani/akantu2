#ifndef __AKANTU_CONTACT_DETECTOR_SHARED_HH__
#define __AKANTU_CONTACT_DETECTOR_SHARED_HH__

#include "aka_common.hh"
#include "aka_grid_dynamic.hh"
#include "fe_engine.hh"

namespace akantu {

/// Base class for contact detectors
class AbstractContactDetector {

  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructors                                                  */
  /* ------------------------------------------------------------------------ */
public:
  AbstractContactDetector(Mesh & mesh, Array<Real> initial_positions);

  virtual ~AbstractContactDetector() = default;

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  // TODO: move body to tmpl header file

  /// Return a new vector with the positions of a node.
  inline Vector<Real> getNodePosition(UInt node) const {
    Vector<Real> position(spatial_dimension);
    for (UInt s : arange(spatial_dimension)) {
      position(s) = this->positions(node, s);
    }
    return position;
  }

  /// Fill the matrix with the coordinates of an element.
  inline void coordinatesOfElement(const Element & el,
                                   Matrix<Real> & coords) const {

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

  /// Compute the minimum and maximum element sizes.
  template <typename vector_type>
  inline auto computeElementSizes(vector_type && nodes) const {
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

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the mesh
  AKANTU_GET_MACRO(Mesh, mesh, Mesh &)

  AKANTU_GET_MACRO_NOT_CONST(Positions, positions, Array<Real> &);
  AKANTU_SET_MACRO(Positions, positions, Array<Real>);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Mesh
  Mesh & mesh;

  /// dimension of the model
  const UInt spatial_dimension;

  /// contains the updated positions of the nodes
  Array<Real> positions; // TODO: get on the fly, storing is wasteful
};

} // namespace akantu

#endif // __AKANTU_CONTACT_DETECTOR_SHARED_HH__
