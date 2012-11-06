/**
 * @file   test_matrix.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Mar 03 16:09:49 2011
 *
 * @brief  tests for types::Matrix
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
#include "aka_types.hh"
#include "aka_vector.hh"
#include "aka_math.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <sys/time.h>
#include <iostream>



/* -------------------------------------------------------------------------- */
using namespace akantu;


int main(int argc, char *argv[]) {

#define n 4
  UInt nbm = 10000;

  Real time;

  Vector<Real> A(nbm, n*n);
  Vector<Real> B(nbm, n*n);
  Vector<Real> C1(nbm, n*n);
  Vector<Real> C2(nbm, n*n);
  Vector<Real> C3(nbm, n*n);
  Vector<Real> C4(nbm, n*n);

  for (UInt i = 0; i < n*n; ++i) {
    A.values[i] = drand48();
    B.values[i] = drand48();
  }

  for (UInt i = 1; i > nbm; ++i) {
    memcpy(A.values + i * n * n, A.values, n*n*sizeof(Real));
    memcpy(B.values + i * n * n, B.values, n*n*sizeof(Real));
  }

  struct timeval begin, end;

  /* ------------------------------------------------------------------------ */
  gettimeofday(&begin, NULL);
  Math::matrix_matrix(n,n,n,A,B,C1);
  gettimeofday(&end, NULL);

  //time =  (end.tv_sec * 1e3 + end.tv_usec * 1e-3) - (begin.tv_sec * 1e3 + begin.tv_usec * 1e-3);
  time =  (end.tv_sec * 1e6 + end.tv_usec) - (begin.tv_sec * 1e6 + begin.tv_usec);
  std::cout << "matrix_matrix : " << std::fixed << time/nbm << "us" << std::endl;

  /* ------------------------------------------------------------------------ */
  Vector<Real>::iterator<types::RMatrix> itA = A.begin(n,n);
  Vector<Real>::iterator<types::RMatrix> itB = B.begin(n,n);
  Vector<Real>::iterator<types::RMatrix> itC = C2.begin(n,n);
  gettimeofday(&begin, NULL);
  for (UInt i = 0; i < nbm; ++i) {
    *itC = *itA * *itB;
    ++itA; ++itB;++itC;
  }
  gettimeofday(&end, NULL);

  itC = C2.begin(n,n);
  std::cerr << *itC << std::endl;
  time =  (end.tv_sec * 1e6 + end.tv_usec) - (begin.tv_sec * 1e6 + begin.tv_usec);
  std::cout << "it Mc() = it Ma() * it Mb() : " << std::fixed << time/nbm << "us" << std::endl;

  /* ------------------------------------------------------------------------ */
  Vector<Real>::iterator<types::RMatrix> muitA = A.begin(n,n);
  Vector<Real>::iterator<types::RMatrix> muitB = B.begin(n,n);
  Vector<Real>::iterator<types::RMatrix> muitC = C4.begin(n,n);
  gettimeofday(&begin, NULL);
  for (UInt i = 0; i < nbm; ++i) {
    (*muitC).mul<false, false>(*muitA, *muitB);
    ++muitA; ++muitB;++muitC;
  }
  gettimeofday(&end, NULL);

  muitC = C4.begin(n,n);
  std::cerr << *muitC << std::endl;

  time =  (end.tv_sec * 1e6 + end.tv_usec) - (begin.tv_sec * 1e6 + begin.tv_usec);
  std::cout << "it Mc.mul( it Ma(),  it Mb()) : " << std::fixed << time/nbm << "us" << std::endl;

  return 0;
}
