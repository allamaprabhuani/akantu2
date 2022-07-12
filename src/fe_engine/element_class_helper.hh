#ifndef ELEMENT_CLASS_HELPER_HH
#define ELEMENT_CLASS_HELPER_HH

#include "element_class.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <ElementKind kind> class ElementClassHelper {};

/* -------------------------------------------------------------------------- */
template <> class ElementClassHelper<_ek_regular> {
public:
  static inline Vector<Real> getN(const Vector<Real> & natural_coords,
                                  ElementType type) {
    return tuple_dispatch<ElementTypes_t<_ek_regular>>(
        [&](auto && enum_type) {
          constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
          constexpr auto nb_nodes_per_element =
              ElementClass<type>::getNbNodesPerElement();
          Vector<Real> shapes(nb_nodes_per_element);
          ElementClass<type>::computeShapes(natural_coords, shapes);
          return shapes;
        },
        type);
  }

  /* ------------------------------------------------------------------------ */
  static inline Matrix<Real> getDNDS(const Vector<Real> & natural_coords,
                                     ElementType type) {
    return tuple_dispatch<ElementTypes_t<_ek_regular>>(
        [&](auto && enum_type) {
          constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
          constexpr auto nb_nodes_per_element =
              ElementClass<type>::getNbNodesPerElement();
          Matrix<Real> dnds(natural_coords.size(), nb_nodes_per_element);
          ElementClass<type>::computeDNDS(natural_coords, dnds);
          return dnds;
        },
        type);
  }

  /* ------------------------------------------------------------------------ */
  static inline Matrix<Real> getD2NDS2(const Vector<Real> & natural_coords,
                                       ElementType type) {
    return tuple_dispatch<ElementTypes_t<_ek_regular>>(
        [&](auto && enum_type) {
          constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
          constexpr auto nb_nodes_per_element =
              ElementClass<type>::getNbNodesPerElement();
          auto dim = natural_coords.size();
          Matrix<Real> d2nds2(dim * dim, nb_nodes_per_element);
          ElementClass<type>::computeD2NDS2(natural_coords, d2nds2);
          return d2nds2;
        },
        type);
  }

  /* ------------------------------------------------------------------------ */
  static inline Matrix<Real> getJMat(const Vector<Real> & natural_coords,
                                     const Matrix<Real> & positions,
                                     ElementType type) {
    return tuple_dispatch<ElementTypes_t<_ek_regular>>(
        [&](auto && enum_type) {
          constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
          constexpr auto nb_nodes_per_element =
              ElementClass<type>::getNbNodesPerElement();
          Matrix<Real> dnds(natural_coords.size(), nb_nodes_per_element);
          Matrix<Real> jmat(dnds.rows(), positions.rows());
          ElementClass<type>::computeDNDS(natural_coords, dnds);
          jmat = ElementClass<type>::computeJMat(dnds, positions);
          return jmat;
        },
        type);
  }
};

} // namespace akantu
#endif // ELEMENT_CLASS_HELPER_H
