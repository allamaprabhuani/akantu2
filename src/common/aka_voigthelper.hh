/**
 * @file   aka_voigthelper
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date   Thu 11 21 13:48:26 2013
 *
 * @brief  Helper file for Voigt notation
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

#ifndef __AKA_VOIGTHELPER_HH__
#define __AKA_VOIGTHELPER_HH__

#include "aka_common.hh"
#include "aka_types.hh"

__BEGIN_AKANTU__

template <UInt dim>
class VoigtHelper {
  public:
    /// transfer the B matrix to a Voigt notation B matrix
    static void transferBMatrixToSymVoigtBMatrix(const Matrix<Real> & B,
          Matrix<Real> & Bvoigt,
          UInt nb_nodes_per_element);

    /// transfer the BNL and BL2 matrix to a Voigt notation B matrix (See Bathe et al. IJNME vol 9, 1975)
    static void transferBMatrixToBNL(const Matrix<Real> & B,
          Matrix<Real> & Bvoigt,
          UInt nb_nodes_per_element);

    static void transferBMatrixToBL2(const Matrix<Real> & B, const Matrix<Real> & grad_u,
          Matrix<Real> & Bvoigt,
          UInt nb_nodes_per_element);

  public:
    const static UInt size;
    // matrix of vector index I as function of tensor indices i,j
    const static UInt mat[dim][dim];
    // array of matrix indices ij as function of vector index I
    const static UInt vec[dim*dim][2];
    // factors to multiply the strain by for voigt notation
    const static Real factors[dim*(dim-(dim-1)/2)];

};

template <UInt dim> const UInt VoigtHelper<dim>::size = dim*(dim-(dim-1)/2);

/* -------------------------------------------------------------------------- */
template<UInt dim>
void VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(const Matrix<Real> & B,
							Matrix<Real> & Bvoigt,
							UInt nb_nodes_per_element) {
  Bvoigt.clear();

  for (UInt i = 0; i < dim; ++i)
    for (UInt n = 0; n < nb_nodes_per_element; ++n)
      Bvoigt(i, i + n*dim) = B(i, n);

  if(dim == 2) {
    ///in 2D, fill the @f$ [\frac{\partial N_i}{\partial x}, \frac{\partial N_i}{\partial y}]@f$ row
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      Bvoigt(2, 1 + n*2) = B(0, n);
      Bvoigt(2, 0 + n*2) = B(1, n);
    }
  }


  if(dim == 3) {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      Real dndx = B(0, n);
      Real dndy = B(1, n);
      Real dndz = B(2, n);

      ///in 3D, fill the @f$ [0, \frac{\partial N_i}{\partial y}, \frac{N_i}{\partial z}]@f$ row
      Bvoigt(3, 1 + n*3) = dndz;
      Bvoigt(3, 2 + n*3) = dndy;

      ///in 3D, fill the @f$ [\frac{\partial N_i}{\partial x}, 0, \frac{N_i}{\partial z}]@f$ row
      Bvoigt(4, 0 + n*3) = dndz;
      Bvoigt(4, 2 + n*3) = dndx;

      ///in 3D, fill the @f$ [\frac{\partial N_i}{\partial x}, \frac{N_i}{\partial y}, 0]@f$ row
      Bvoigt(5, 0 + n*3) = dndy;
      Bvoigt(5, 1 + n*3) = dndx;
    }
  }
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void VoigtHelper<dim>::transferBMatrixToBNL(const Matrix<Real> & B,
        Matrix<Real> & Bvoigt,
        UInt nb_nodes_per_element) {

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

__END_AKANTU__

#endif
