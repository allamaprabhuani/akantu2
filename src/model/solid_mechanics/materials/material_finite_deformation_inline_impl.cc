/**
 * @file   material_finite_deformation_inline_impl.cc
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Wed Aug 04 10:58:42 2010
 *
 * @brief  Implementation of the inline functions of the material elastic
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

#include <iostream>

#include "aka_error.hh"


/* -------------------------------------------------------------------------- */

/*************************************************************************************************************/

/*template<UInt dim>
inline void MaterialFiniteDeformation::transferBMatrixToBL2(const Matrix<Real> & B, const Matrix<Real> & grad_u,
        Matrix<Real> & Bvoigt,
        UInt nb_nodes_per_element) const {

    Bvoigt.clear();

    //see Finite ekement formulations for large deformation dynamic analysis, Bathe et al. IJNME vol 9, 1975, page 364 B_{NL}

    for (UInt i = 0; i < 3 * (dim - 1); ++i) {
        for (UInt j = 0; j < nb_nodes_per_element; ++j) {
            for (UInt k = 0; k < dim; ++k) {
                if (i < dim)
                  Bvoigt(i, j * dim + k) = grad_u(k, i) * B(i, j);
                else {
                    if (dim == 3) {
                        UInt aux = i ;
                        if (aux >= dim)
                            aux -= dim;
                        for (UInt m = 0; m < dim; ++m) {
                            if (m != aux) {
                                UInt index1 = m;
                                UInt index2 = dim - m - aux;
                                Bvoigt(i, j * dim + k) += grad_u(k, index1) * B(index2, j);
                            }
                        }
                    } else if (dim == 2) {
                        for (UInt m = 0; m < dim; ++m) {
                            UInt index1 = m;
                            UInt index2 = (dim - 1) - m;
                            Bvoigt(i, j * dim + k) += grad_u(k, index1) * B(index2, j);
                        }
                    }

                }
            }
        }

        //TODO: Verify the 2D and 1D case
    }
    }*/

/*template<UInt dim>
inline void MaterialFiniteDeformation::transferBMatrixToBL1(const Matrix<Real> & B,
        Matrix<Real> & Bvoigt,
        UInt nb_nodes_per_element) const {

    Bvoigt.clear();

    //see Finite ekement formulations for large deformation dynamic analysis, Bathe et al. IJNME vol 9, 1975, page 364 B_{NL}

    for (UInt i = 0; i < 3 * (dim - 1); ++i) {
        for (UInt j = 0; j < nb_nodes_per_element; ++j) {
            if (i < dim)
              Bvoigt(i, j * dim + i) = B(i, j);
            else {
                if (dim == 3) {
                    UInt aux = i - 1;
                    if (aux >= dim)
                        aux -= dim;
                    for (UInt k = 0; k < dim; ++k) {
                        if (k != aux) {
                            UInt index = dim - k - aux;
                            Bvoigt(i, j * dim + k) = B(index, j);
                        }
                    }
                } else if (dim == 2) {
                    for (UInt k = 0; k < dim; ++k) {
                        UInt index = (dim - 1) - k;
                        Bvoigt(i, j * dim + k) = B(index, j);

                    }

                }
            }
        }
    }

    //TODO: Verify the 2D and 1D case
}*/

/* -------------------------------------------------------------------------- */
/*template<UInt dim>
inline void MaterialFiniteDeformation::GreenStrain(const Matrix<Real> & F,
        Matrix<Real> & E) {

    E.mul < true, false > (F, F);

    for (UInt i = 0; i < dim; ++i) {
        E(i, i) -= 1.0;
        for (UInt j = 0; j < dim; ++j)
            E(i, j) *= 0.5;
    }
}*/

/* -------------------------------------------------------------------------- */
/*template<UInt dim>
inline void MaterialFiniteDeformation::deltaGreenStrain(const Matrix<Real> & d_u, const Matrix<Real> & u,
        Matrix<Real> & E) {

    Matrix<Real> F_tensor(3, 3);
    Matrix<Real> E_tensor(3, 3);
    Matrix<Real> Aux(3, 3);
    Matrix<Real> Aux2(3, 3);

    gradUToF<dim > (d_u, F_tensor);
    GreenStrain<dim > (F_tensor, E_tensor);

    Aux.mul < true, false > (d_u, u);
    Aux2.mul < true, false > (u, d_u);

    for (UInt i = 0; i < dim; ++i)
        for (UInt j = 0; j < dim; ++j)
            E(i, j) = E_tensor(i, j) + 0.5 * (Aux(i, j) + Aux2(i, j));

            }*/
