/**
 * @file   element_class_line_1.cc
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
template<> inline void ElementClass<_line_1>::shapeFunctions(const Real * x,
							     Real * shape,
							     Real * shape_deriv,
							     Real * jacobian) {
  /**
   *           q
   *  x--------|--------x
   * -1        0        1
   *
   * N1(e) = 1/2 (1 - e)
   * N2(e) = 1/2 (1 + e)
   */

  Real weight = 2;

  /// shape functions
  shape[0] = .5; //N1(0)
  shape[1] = .5; //N2(0)

  Real dnds[nb_nodes_per_element*spatial_dimension];
  ///dN1/de
  dnds[0] = - .5;
  ///dN2/de
  dnds[1] =   .5;

  Real dxds = dnds[0]*x[0] + dnds[1]*x[1];

  jacobian[0] = dxds * weight;

  shape_deriv[0] = dnds[0] / dxds;
  shape_deriv[1] = dnds[1] / dxds;
}


/* ------------------------------------------------------------------------ */
template<> inline Real ElementClass<_line_1>::getVolume(const double * coord) {
  return sqrt((coord[0] - coord[1])*(coord[0] - coord[1])); //*1*1 for volume
}

/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_line_1>::getInradius(const double * coord) {
  return sqrt((coord[0] - coord[1])*(coord[0] - coord[1]));
}
