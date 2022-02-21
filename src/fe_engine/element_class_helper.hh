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
#define GET_SHAPE_NATURAL(type)                                                \
  auto nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();      \
  Vector<Real> shapes(nb_nodes_per_element);                                   \
  ElementClass<type>::computeShapes(natural_coords, shapes);                   \
  return shapes
    AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_SHAPE_NATURAL);
#undef GET_SHAPE_NATURAL

    return Vector<Real>(0);
  }

  /* ------------------------------------------------------------------------ */
  static inline Matrix<Real> getDNDS(const Vector<Real> & natural_coords,
                                     ElementType type) {
#define GET_DNDS_NATURAL(type)                                                 \
  auto nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();      \
  Matrix<Real> dnds(natural_coords.size(), nb_nodes_per_element);              \
  ElementClass<type>::computeDNDS(natural_coords, dnds);                       \
  return dnds
    AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_DNDS_NATURAL);
#undef GET_DNDS_NATURAL

    return Matrix<Real>(0, 0);
  }

  /* ------------------------------------------------------------------------ */
  static inline Matrix<Real> getD2NDS2(const Vector<Real> & natural_coords,
                                       ElementType type) {
#define GET_D2ND2S_NATURAL(type)                                               \
  auto nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();      \
  auto dim = natural_coords.size();                                            \
  Matrix<Real> d2nds2(dim * dim, nb_nodes_per_element);                        \
  ElementClass<type>::computeD2NDS2(natural_coords, d2nds2);                   \
  return d2nds2
    AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_D2ND2S_NATURAL);
#undef GET_D2NDS2_NATURAL

    return Matrix<Real>(0, 0);
  }

  /* ------------------------------------------------------------------------ */
  static inline Matrix<Real> getJMat(const Vector<Real> & natural_coords,
                                     const Matrix<Real> & positions,
                                     ElementType type) {
#define GET_JMAT_NATURAL(type)                                                 \
  auto nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement();      \
  Matrix<Real> dnds(natural_coords.size(), nb_nodes_per_element);              \
  Matrix<Real> jmat(dnds.rows(), positions.rows());                            \
  ElementClass<type>::computeDNDS(natural_coords, dnds);                       \
  ElementClass<type>::computeJMat(dnds, positions, jmat);                      \
  return jmat
    AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_JMAT_NATURAL);
#undef GET_JMAT_NATURAL

    return Matrix<Real>(0, 0);
  }
};

} // namespace akantu
#endif // ELEMENT_CLASS_HELPER_H
