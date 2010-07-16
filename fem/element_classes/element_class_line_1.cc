/**
 * @file   element_class_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 10:28:28 2010
 *
 * @brief  Specialization of the element_class class for the type _line_1
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
template<> ElementClass<_line_1>::ElementClass() {
  nb_nodes_per_element = 2;
  nb_quadrature_points = 1;
  spatial_dimension    = 1;
}

/* -------------------------------------------------------------------------- */
template<> void ElementClass<_line_1>::shapeFunctions(const Real * x,
						      Real * shape,
						      Real * dshape,
						      Real * jacobian) {
  /**
   *           q
   *  x--------|--------x
   * -1        0        1
   *
   * N1(e) = 1/2 (1 - e)
   * N2(e) = 1/2 (1 + e)
   */

  /// shape functions
  shape[0] = .5;
  shape[1] = .5;


  Real dnds[nb_nodes_per_element*spatial_dimension];
  ///dN1/de
  dnds[0] = - .5;
  ///dN2/de
  dnds[1] =   .5;

  Real dxds = dnds[0]*x[0] + dnds[1]*x[1];

  jacobian[0] = 1 / dxds;

  dshape[0] = dnds[0] / dxds;
  dshape[1] = dnds[1] / dxds;

}


/* -------------------------------------------------------------------------- */
template<> Real ElementClass<_line_1>::volume(const double * coord) {
  return sqrt((coord[0] - coord[1])*(coord[0] - coord[1])); // * 1 * 1 for volume;
}

/* -------------------------------------------------------------------------- */
