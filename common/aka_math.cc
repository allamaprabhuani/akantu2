/**
 * @file   aka_math.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jul 28 12:13:46 2010
 *
 * @brief  Implementation of the math toolbox
 *
 * @section LICENSE
 *
 * \<insert license here\>
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
			 Vector<Real> & y) {
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
    matrix_vector(m, n, A_val, x_val, y_val);

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
			 Vector<Real> & C) {
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
    matrix_matrix(m, n, k, A_val, B_val, C_val);

    A_val += offset_A;
    B_val += offset_B;
    C_val += offset_C;
  }

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
