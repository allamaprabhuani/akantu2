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
inline UInt MaterialFiniteDeformation::getCauchyStressMatrixSize(UInt dim) const {
    //switch (dim){ case 1: return 1; case 2: return 4; case 3: return 9;}
    //return (dim * dim + (3 - dim));
    return (dim * dim);
}

/* -------------------------------------------------------------------------- */
inline UInt MaterialFiniteDeformation::getCauchyStressArraySize(UInt dim) const {
    switch(dim){case 1: return 3; case 2: return 3; case 3: return 6; default: AKANTU_DEBUG_ERROR("The dimension is not right!!");};
    return 0;
}

/* -------------------------------------------------------------------------- */


template<UInt dim>
inline void MaterialFiniteDeformation::transferBMatrixToBNL(const Matrix<Real> & B,
        Matrix<Real> & Bvoigt,
        UInt nb_nodes_per_element) const {

    Bvoigt.clear();

    //see Finite element formulations for large deformation dynamic analysis, Bathe et al. IJNME vol 9, 1975, page 364 B_{NL}

    for (UInt i = 0; i < dim; ++i) {
        for (UInt m = 0; m < nb_nodes_per_element; ++m) {
            for (UInt n = 0; n < dim; ++n) {
                //std::cout << B(n, m) << std::endl;
	      Bvoigt(i * dim + n, m * dim + i) = B(n, m);
            }
        }
    }

    //TODO: Verify the 2D and 1D case
}

template<UInt dim>
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
                        UInt aux = i - 1;
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
}

template<UInt dim>
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
}

template<UInt dim>
inline void MaterialFiniteDeformation::SetCauchyStressMatrix(const Matrix<Real> & S_t, Matrix<Real> & Stress_matrix) {

    AKANTU_DEBUG_IN();

    Stress_matrix.clear();

    //see Finite ekement formulations for large deformation dynamic analysis, Bathe et al. IJNME vol 9, 1975, page 364 ^t\tau

    for (UInt i = 0; i < dim; ++i) {
        for (UInt m = 0; m < dim; ++m) {
            for (UInt n = 0; n < dim; ++n) {
                Stress_matrix(i * dim + m, i * dim + n) = S_t(m, n);
            }
        }
    }

    //other terms from the diagonal
    /*for (UInt i = 0; i < 3 - dim; ++i) {
        Stress_matrix(dim * dim + i, dim * dim + i) = S_t(dim + i, dim + i);
    }*/


    AKANTU_DEBUG_OUT();
}

/* ---------------------------------------------------------------------------*/
template<UInt dim>
inline void MaterialFiniteDeformation::SetCauchyStressArray(const Matrix<Real> & S_t, Matrix<Real> & Stress_vect) {

    AKANTU_DEBUG_IN();

    Stress_vect.clear();
    
    //UInt cauchy_matrix_size = getCauchyStressArraySize(dim);

    //see Finite ekement formulations for large deformation dynamic analysis, Bathe et al. IJNME vol 9, 1975, page 364 ^t\tau

    /*
     * 1d: [ s11 ]'
     * 2d: [ s11 s22 s12 ]'
     * 3d: [ s11 s22 s33 s12 s23 s13 ]
     */
    for (UInt i = 0; i < dim; ++i)//diagonal terms
        Stress_vect(i, 0) = S_t(i, i);
    
    for (UInt i = 1; i < dim; ++i)// term s12 in 2D and terms s12 s23 in 3D
        Stress_vect(dim+i-1, 0) = S_t(i-1, i);
    
    for (UInt i = 2; i < dim; ++i)//term s13 in 3D
        Stress_vect(dim+i, 0) = S_t(i-2, i);
    

    /* Append some values in 1d and 2d:
     * 1d: [ s22 s33 ]'
     * 2d: [ s33 ]
     */
    /*for (UInt i = 0; i < 3 - dim; ++i) 
        Stress_vect(dim * dim + i, 0) = S_t(dim + i, dim + i);*/

    AKANTU_DEBUG_OUT();
}

/* ---------------------------------------------------------------------------*/
template<UInt dim>
inline void MaterialFiniteDeformation::SetCauchyStrainArray(const Matrix<Real> & E_t, Matrix<Real> & Strain_vect) {

    AKANTU_DEBUG_IN();

    Strain_vect.clear();
    
    //UInt cauchy_matrix_size = getCauchyStressArraySize(dim);

    //see Finite ekement formulations for large deformation dynamic analysis, Bathe et al. IJNME vol 9, 1975, page 364 ^t\tau

    /*
     * 1d: [ s11 ]'
     * 2d: [ s11 s22 s12 ]'
     * 3d: [ s11 s22 s33 s12 s23 s13 ]
     */
    for (UInt i = 0; i < dim; ++i)//diagonal terms
        Strain_vect(i, 0) = E_t(i, i);
    
    for (UInt i = 1; i < dim; ++i)// term s12 in 2D and terms s12 s23 in 3D
        Strain_vect(dim+i-1, 0) = E_t(i-1, i)+E_t(i, i-1);
    
    for (UInt i = 2; i < dim; ++i)//term s13 in 3D
        Strain_vect(dim+i, 0) = E_t(i-2, i) + E_t(i, i-2);
    

    /* Append some values in 1d and 2d:
     * 1d: [ s22 s33 ]'
     * 2d: [ s33 ]
     */
    for (UInt i = 0; i < 3 - dim; ++i) 
        Strain_vect(dim * dim + i, 0) = E_t(dim + i, dim + i);

    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void MaterialFiniteDeformation::GreenStrain(const Matrix<Real> & F,
        Matrix<Real> & E) {

    E.mul < true, false > (F, F);

    for (UInt i = 0; i < dim; ++i) {
        E(i, i) -= 1.0;
        for (UInt j = 0; j < dim; ++j)
            E(i, j) *= 0.5;
    }
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void MaterialFiniteDeformation::AlmansiStrain(const Matrix<Real> & F,
        Matrix<Real> & E) {
    
    Matrix<Real> Fminus(3, 3);
    Math::inv3(F.storage(),Fminus.storage());
    E.mul < true, false > (Fminus, Fminus);

    for (UInt i = 0; i < dim; ++i) {
        E(i, i) -= 1.0;
        for (UInt j = 0; j < dim; ++j)
            E(i, j) *= -0.5;
    }
    
}
/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void MaterialFiniteDeformation::GreenStrain(const Matrix<Real> & d_u, const Matrix<Real> & u,
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

}
