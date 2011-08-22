/**
 * @file   element_class_bernoulli_beam_2.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu Mar 31 14:02:22 2011
 *
 * @brief  Specialization of the element_class class for the type _bernoulli_beam_2
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 * @section DESCRIPTION
 *
 * @verbatim

   --x-----q1----|----q2-----x---> x
    -a          0            a
 @endverbatim
 *
 * @subsection coords Nodes coordinates
 *
 * @f[
 * \begin{array}{lll}
 *  x_{1}  = -a   &  x_{2} = a
 * \end{array}
 * @f]
 *
 * @subsection shapes Shape functions

 * @f[
 * \begin{description}
 * \item[id=0]
 *
 * \begin{array}{ll}
 * N_1(x) = 1/2(1-x/a)\\
 * N_2(x) = 1/2(1+x/a)\\
 * \end{array}
 *
 * \item[id=1]
 *
 * \begin{array}{ll}
 * M_1(x) = 1/4(x^{3}/a^{3}-3x/a+2)\\
 * M_2(x) = -1/4(x^{3}/a^{3}-3x/a-2)\\
 * \end{array}
 *
 * \item[id=2]
 *
 * \begin{array}{ll}
 * L_1(x) = a/4(x^{3}/a^{3}-x^{2}/a^{2}-x/a+1)\\
 * L_2(x) = a/4(x^{3}/a^{3}+x^{2}/a^{2}-x/a-1)\\
 * \end{array}
 *
 * \item[id=3]
 *
 * \begin{array}{ll}
 * M'_1(x) = 3/4a(x^{2}/a^{2}-1)\\
 * M'_2(x) = -3/4a(x^{2}/a^{2}-1)\\
 * \end{array}
 *
 *\item[id=4]
 *
 *\begin{array}{ll}
 * L'_1(x) = 1/4(3x^{2}/a^{2}-2x/a-1)\\
 * L'_2(x) = 1/4(3x^{2}/a^{2}+2x/a-1)\\
 *\end{array}
 *\end{description}
 *@f]
 *
 * @subsection dnds Shape derivatives
 *
 *@f[
 *\begin{description}
 *\item[id=0]
 *\begin{array}{ll}
 *N'_1(x) = -1/2a\\
 *N'_2(x) = 1/2a\\
 *\end{array}]
 *\item[id=1]
 *\begin{array}{ll}
 *-M''_1(x) = -3x/(2a^{3})\\
 *-M''_2(x) = 3x/(2a^{3})\\
 *\end{array}
 *\item[id=2]
 *\begin{array}{ll}
 *-L''_1(x) = -1/2a(3x/a-1)\\
 *-L''_2(x) = -1/2a(3x/a+1)\\
 *\end{array}
 *\end{description}
 *@f]
 *
 * @subsection quad_points Position of quadrature points
 *
 * @f[
 * \begin{array}{ll}
 * x_{q1}  = -a/\sqrt{3} & x_{q2} = a/\sqrt{3}
 * \end{array}
 * @f]
 */

/* -------------------------------------------------------------------------- */


template<> UInt ElementClass<_bernoulli_beam_2>::nb_nodes_per_element;
template<> UInt ElementClass<_bernoulli_beam_2>::nb_quadrature_points;
template<> UInt ElementClass<_bernoulli_beam_2>::spatial_dimension;

/* -------------------------------------------------------------------------- */
template <>
inline void ElementClass<_bernoulli_beam_2>::computeShapes(const Real * natural_coords,
							   Real * shapes,
							   const Real * local_coord,
							   UInt id) {
  /// Compute the dimension of the beam
  Real a=.5 * Math::distance_2d(local_coord, local_coord+2);

  /// natural coordinate
  Real c =(*natural_coords)*a;


  switch (id) {

  case 0:
    shapes[0]=0.5*(1-c/a);
    shapes[1]=0.5*(1+c/a);
    break;

  case 1:
    shapes[0]=0.25*(pow(c,3)/pow(a,3)-3*c/a+2);
    shapes[1]=-0.25*(pow(c,3)/pow(a,3)-3*c/a-2);
    break;

  case 2:
    shapes[0]=0.25*a*(pow(c,3)/pow(a,3)-c*c/(a*a)-c/a+1);
    shapes[1]=0.25*a*(pow(c,3)/pow(a,3)+c*c/(a*a)-c/a-1);
    break;

  case 3:
    shapes[0]=0.75/a*(c*c/(a*a)-1);
    shapes[1]=-0.75/a*(c*c/(a*a)-1);
  break;

  case 4:
    shapes[0]=0.25*(3*c*c/(a*a)-2*c/a-1);
    shapes[1]=0.25*(3*c*c/(a*a)+2*c/a-1);
  break;
  }
 }

/* -------------------------------------------------------------------------- */
template <>
inline void ElementClass<_bernoulli_beam_2>::computeShapeDerivatives(const Real * natural_coords,
								     Real * shape_deriv,
								     const Real * local_coord,
								     UInt id) {


 /// Compute the dimension of the beam
  Real a=0.5*Math::distance_2d(local_coord,local_coord+2);
  Real x1=*local_coord;
  Real y1=*(local_coord+1);
  Real x2=*(local_coord+2);
  Real y2=*(local_coord+3);
  Real tetha=std::atan((y2-y1)/(x2-x1));
  Real pi = std::atan(1.0)*4;

  if ((x2-x1) < 0) {
    tetha += pi;
  }

  /// natural coordinate
  Real c = (*natural_coords)*a;

  /// Definition of the rotation matrix
  Real T[4];
  T[0]= cos(tetha);
  T[1]= sin(tetha);
  T[2]=-T[1];
  T[3]= T[0];

  // B archetype
  Real shape_deriv_arch[4];
  switch (id) {

  case 0:
    shape_deriv_arch[0]=-0.5/a;
    shape_deriv_arch[1]= 0;
    shape_deriv_arch[2]= 0.5/a;
    shape_deriv_arch[3]= 0;

    Math::matrix_matrix(2,2,2,shape_deriv_arch,T,shape_deriv);

    break;

  case 1:
    shape_deriv_arch[0]= 0.;
    shape_deriv_arch[1]=-3.*c/(2.*pow(a,3));
    shape_deriv_arch[2]= 0.;
    shape_deriv_arch[3]= 3.*c/(2.*pow(a,3));

    Math::matrix_matrix(2,2,2,shape_deriv_arch,T,shape_deriv);

    break;

  case 2:
    shape_deriv[0]=-0.5/a*(3*c/a-1);
    shape_deriv[1]=-0.5/a*(3*c/a+1);
    shape_deriv[2]=0;
    shape_deriv[3]=0;
    break;
  }

}

/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_bernoulli_beam_2>::computeJacobian(const Real * coord,
									 const UInt nb_points,
									 __attribute__((unused)) const UInt dimension,
									 Real * jac){
  Real a=0.5*Math::distance_2d(coord,coord+2);

  for (UInt p = 0; p < nb_points; ++p) {
    jac[p]=a;
  }
}

/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_bernoulli_beam_2>::getInradius(const Real * coord) {
   return Math::distance_2d(coord, coord+2);
}

/* -------------------------------------------------------------------------- */
