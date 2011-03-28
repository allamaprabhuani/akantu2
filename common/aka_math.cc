/**
 * @file   aka_math.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jul 28 12:13:46 2010
 *
 * @brief  Implementation of the math toolbox
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
 */

/* -------------------------------------------------------------------------- */
#include "aka_math.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void Math::matrix_vector(UInt m, UInt n,
			 const Vector<Real> & A,
			 const Vector<Real> & x,
			 Vector<Real> & y,
			 Real alpha) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(A.getSize() == x.getSize(),
		      "The vector A(" << A.getID()
		      << ") and the vector x(" << x.getID()
		      << ") must have the same size");

  AKANTU_DEBUG_ASSERT(A.getNbComponent() == m * n,
		      "The vector A(" << A.getID()
		      << ") has the good number of component.");

  AKANTU_DEBUG_ASSERT(x.getNbComponent() == n,
		      "The vector x(" << x.getID()
		      << ") do not the good number of component.");

  AKANTU_DEBUG_ASSERT(y.getNbComponent() == n,
		      "The vector y(" << y.getID()
		      << ") do not the good number of component.");


  UInt nb_element = A.getSize();
  UInt offset_A = A.getNbComponent();
  UInt offset_x = x.getNbComponent();

  y.resize(nb_element);

  Real * A_val = A.values;
  Real * x_val = x.values;
  Real * y_val = y.values;

  for (UInt el = 0; el < nb_element; ++el) {
    matrix_vector(m, n, A_val, x_val, y_val, alpha);

    A_val += offset_A;
    x_val += offset_x;
    y_val += offset_x;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Math::matrix_matrix(UInt m, UInt n, UInt k,
			 const Vector<Real> & A,
			 const Vector<Real> & B,
			 Vector<Real> & C,
			 Real alpha) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(A.getSize() == B.getSize(),
		      "The vector A(" << A.getID()
		      << ") and the vector B(" << B.getID()
		      << ") must have the same size");

  AKANTU_DEBUG_ASSERT(A.getNbComponent() == m * k,
		      "The vector A(" << A.getID()
		      << ") has the good number of component.");

  AKANTU_DEBUG_ASSERT(B.getNbComponent() == k * n ,
		      "The vector B(" << B.getID()
		      << ") do not the good number of component.");

  AKANTU_DEBUG_ASSERT(C.getNbComponent() == m * n,
		      "The vector C(" << C.getID()
		      << ") do not the good number of component.");

  UInt nb_element = A.getSize();
  UInt offset_A = A.getNbComponent();
  UInt offset_B = B.getNbComponent();
  UInt offset_C = C.getNbComponent();

  C.resize(nb_element);

  Real * A_val = A.values;
  Real * B_val = B.values;
  Real * C_val = C.values;

  for (UInt el = 0; el < nb_element; ++el) {
    matrix_matrix(m, n, k, A_val, B_val, C_val, alpha);

    A_val += offset_A;
    B_val += offset_B;
    C_val += offset_C;
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void Math::matrix_matrixt(UInt m, UInt n, UInt k,
			  const Vector<Real> & A,
			  const Vector<Real> & B,
			  Vector<Real> & C,
			  Real alpha) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(A.getSize() == B.getSize(),
		      "The vector A(" << A.getID()
		      << ") and the vector B(" << B.getID()
		      << ") must have the same size");

  AKANTU_DEBUG_ASSERT(A.getNbComponent() == m * k,
		      "The vector A(" << A.getID()
		      << ") has the good number of component.");

  AKANTU_DEBUG_ASSERT(B.getNbComponent() == k * n ,
		      "The vector B(" << B.getID()
		      << ") do not the good number of component.");

  AKANTU_DEBUG_ASSERT(C.getNbComponent() == m * n,
		      "The vector C(" << C.getID()
		      << ") do not the good number of component.");

  UInt nb_element = A.getSize();
  UInt offset_A = A.getNbComponent();
  UInt offset_B = B.getNbComponent();
  UInt offset_C = C.getNbComponent();

  C.resize(nb_element);

  Real * A_val = A.values;
  Real * B_val = B.values;
  Real * C_val = C.values;

  for (UInt el = 0; el < nb_element; ++el) {
    matrix_matrixt(m, n, k, A_val, B_val, C_val, alpha);

    A_val += offset_A;
    B_val += offset_B;
    C_val += offset_C;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Math::matrix33_eigenvalues(Real * A, Real *Adiag) {
  ///d = determinant of Matrix A
  Real d = det3(A);
  ///b = trace of Matrix A
  Real b = A[0]+A[4]+A[8];

  Real a = -1 ;
  /// c = 0.5*(trace(M^2)-trace(M)^2)
  Real c =  A[3]*A[1] + A[2]*A[6] + A[5]*A[7] - A[0]*A[4] -
    A[0]*A[8] - A[4]*A[8];
  /// Define x, y, z
  Real x = ((3*c/a) - ((b*b)/(a*a)))/3;
  Real y=((2*(b*b*b)/(a*a*a)) - (9*b*c/(a*a)) + (27*d/a))/27;
  Real z = (y*y)/4 + (x*x*x)/27;
  /// Define I, j, k, m, n, p (so equations are not so cluttered)
  Real i = sqrt(y*y/4 - z);
  Real j = -pow(i,1./3.);
  Real k;
  if (fabs(i)<1e-12)
     k = 0;
     else
     k = acos(-(y/(2*i)));

  Real m = cos(k/3);
  Real n = sqrt(3)*sin(k/3);
  Real p = b/(3*a);

  Adiag[0]=-(2*j*m + p);;
  Adiag[1]=-(-j *(m + n) + p);
  Adiag[2]=-(-j * (m - n) + p);
}


__END_AKANTU__
