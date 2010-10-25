/**
 * @file   element_class_line_2.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Sun Oct 3 10:28:28 2010
 *
 * @brief  Specialization of the element_class class for the type _line_2
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 * @section DESCRIPTION
 *
 * @verbatim

             -1         0         1
  	 -----x---------x---------x-----> x
	      1         3         2
 @endverbatim
 *
 * @subsection coords Nodes coordinates
 *
 * @f[
 * \begin{array}{lll}
 *   x_{1}  = -1   &  x_{2} = 1 & x_{3} = 0
 * \end{array}
 * @f]
 *
 * @subsection shapes Shape functions
 * @f[
 * \begin{array}{ll}
 *   w_1(x) = \frac{x}{2}(x - 1) & w'_1(x) = x - \frac{1}{2}\\
 *   w_2(x) = \frac{x}{2}(x + 1) & w'_2(x) = x + \frac{1}{2}\\
 *   w_3(x) =  1-x^2 & w'_3(x) = -2x
 * \end{array}
 * @f]
 *
 * @subsection quad_points Position of quadrature points
 * @f[
 * \begin{array}{ll}
 * x_{q1}  = -1/\sqrt{3} & x_{q2} = 1/\sqrt{3}
 * \end{array}
 * @f]
 */

  // /// quadrature point position
  // Real quad[spatial_dimension * nb_quadrature_points];
  // quad[0] = -1./sqrt(3);
  // quad[1] = 1./sqrt(3);


/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_line_2>::nb_nodes_per_element;
template<> UInt ElementClass<_line_2>::nb_quadrature_points;
template<> UInt ElementClass<_line_2>::spatial_dimension;



/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_line_2>::computeShapes(const Real * natural_coords, 
							    Real * shapes){
  Real c = natural_coords[0];
  shapes[0] = (c - 1) * c / 2;
  shapes[1] = (c + 1) * c / 2;
  shapes[2] = 1 - c * c;
}
/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_line_2>::computeDNDS(const Real * natural_coords,
							  Real * dnds){

  Real c = natural_coords[0];
  dnds[0]  = c - .5;
  dnds[1]  = c + .5;
  dnds[2]  = -2 * c;
}


/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_line_2>::computeJacobian(const Real * dxds,
							       const UInt dimension, 
							       Real & jac){
  if (dimension == spatial_dimension){
    Real weight = 1;
    jac = dxds[0] * weight;
  }  
  else {
    jac = Math::norm2(dxds);
  }
}

/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_line_2>::getInradius(const Real * coord) {
  Real dist1 = sqrt((coord[0] - coord[1])*(coord[0] - coord[1]));
  Real dist2 = sqrt((coord[1] - coord[2])*(coord[1] - coord[2]));
  return dist1 < dist2 ? dist1 : dist2;
}
