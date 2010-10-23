/**
 * @file   element_class_line_1.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 10:28:28 2010
 *
 * @brief  Specialization of the element_class class for the type _line_1
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 * @section DESCRIPTION
 *
 * @verbatim
              q
   --x--------|--------x---> x
    -1        0        1
 @endverbatim
 *
 * @subsection shapes Shape functions
 * @f{eqnarray*}{
 * w_1(x) &=& 1/2(1 - x) \\
 * w_2(x) &=& 1/2(1 + x)
 * @f}
 *
 * @subsection quad_points Position of quadrature points
 * @f{eqnarray*}{
 * x_{q}  &=& 0
 * @f}
 */

/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_line_1>::nb_nodes_per_element;
template<> UInt ElementClass<_line_1>::nb_quadrature_points;
template<> UInt ElementClass<_line_1>::spatial_dimension;


/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_line_1>::computeShapes(const Real * natural_coords, 
							     Real * shapes){

  /// natural coordinate
  Real c = natural_coords[0];
  /// shape functions
  shapes[0] = 0.5*(1-c);
  shapes[1] =0.5*(1+c);
}
/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_line_1>::computeDNDS(__attribute__ ((unused)) 
							   const Real * natural_coords,
							   Real * dnds){
  /// dN1/de
  dnds[0] = - .5;
  /// dN2/de
  dnds[1] =   .5;
}


/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_line_1>::computeJacobian(const Real * dxds,
								   const UInt dimension, 
								   Real & jac){

  if (dimension == spatial_dimension){
    Real weight = 2.;
    jac = dxds[0]*weight;
  }
  else {
    AKANTU_DEBUG_ERROR("to be implemented");
  }
}
 
/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_line_1>::getInradius(const Real * coord) {
  return sqrt((coord[0] - coord[1])*(coord[0] - coord[1]));
}

