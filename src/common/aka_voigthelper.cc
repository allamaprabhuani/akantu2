/**
 * @file   aka_voigthelper
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date   Thu 11 21 15:25:59 2013
 *
 * @brief  Voigt indices
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2013 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "aka_common.hh"
#include "aka_voigthelper.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */

template <> const UInt VoigtHelper<1>::mat[][1] = {{0}};
template <> const UInt VoigtHelper<2>::mat[][2] = {{0, 2},
                                                  {3, 1}};
template <> const UInt VoigtHelper<3>::mat[][3] = {{0, 5, 4},
                                                   {8, 1, 3},
                                                   {7, 6, 2}};
template <> const UInt VoigtHelper<1>::vec[][2] = {{0, 0}};
template <> const UInt VoigtHelper<2>::vec[][2] = {{0, 0},
                                                   {1, 1},
                                                   {0, 1},
                                                   {1, 0}};
template <> const UInt VoigtHelper<3>::vec[][2] = {{0, 0},
                                                   {1, 1},
                                                   {2, 2},
                                                   {1, 2},
                                                   {0, 2},
                                                   {0, 1},
                                                   {2, 1},
                                                   {2, 0},
                                                   {1, 0}};
template <> const Real VoigtHelper<1>::factors[] = {1.};
template <> const Real VoigtHelper<2>::factors[] = {1., 1., 1., 2.};
template <> const Real VoigtHelper<3>::factors[] = {1., 1., 1.,
                                                    2., 2., 2.};

/* -------------------------------------------------------------------------- */
template<>
void VoigtHelper<3>::transferBMatrixToBL2(const Matrix<Real> & B, const Matrix<Real> & grad_u,
        Matrix<Real> & Bvoigt,
        UInt nb_nodes_per_element) {

    Bvoigt.clear();

    for (UInt i = 0; i < 3; ++i)
        for (UInt j = 0; j < nb_nodes_per_element; ++j)
            for (UInt k = 0; k < 3; ++k)
                  Bvoigt(i, j * 3 + k) = grad_u(k, i) * B(i, j);

    for (UInt i = 3; i < 6; ++i)
        for (UInt j = 0; j < nb_nodes_per_element; ++j)
          for (UInt k = 0; k < 3; ++k){
            UInt aux = i-3 ;
            for (UInt m = 0; m < 3; ++m) {
              if (m != aux) {
                UInt index1 = m;
                UInt index2 = 3 - m - aux;
                Bvoigt(i, j * 3 + k) += grad_u(k, index1) * B(index2, j);
              }
            }
          }
}

/* -------------------------------------------------------------------------- */
template<>
void VoigtHelper<2>::transferBMatrixToBL2(const Matrix<Real> & B, const Matrix<Real> & grad_u,
        Matrix<Real> & Bvoigt,
        UInt nb_nodes_per_element) {

    Bvoigt.clear();

    for (UInt i = 0; i < 2; ++i)
        for (UInt j = 0; j < nb_nodes_per_element; ++j)
            for (UInt k = 0; k < 2; ++k)
                  Bvoigt(i, j * 2 + k) = grad_u(k, i) * B(i, j);

    for (UInt j = 0; j < nb_nodes_per_element; ++j)
      for (UInt k = 0; k < 2; ++k)
        for (UInt m = 0; m < 2; ++m) {
          UInt index1 = m;
          UInt index2 = (2 - 1) - m;
          Bvoigt(2, j * 2 + k) += grad_u(k, index1) * B(index2, j);
        }
}

/* -------------------------------------------------------------------------- */
template<>
void VoigtHelper<1>::transferBMatrixToBL2(const Matrix<Real> & B, const Matrix<Real> & grad_u,
        Matrix<Real> & Bvoigt,
        UInt nb_nodes_per_element) {

    Bvoigt.clear();

    for (UInt j = 0; j < nb_nodes_per_element; ++j)
      for (UInt k = 0; k < 2; ++k)
        Bvoigt(0, j * 2 + k) = grad_u(k, 0) * B(0, j);
}

__END_AKANTU__

