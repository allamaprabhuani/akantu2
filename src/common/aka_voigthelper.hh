
/**
 * @file   aka_voigthelper.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Till Junge <till.junge@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Dec 20 2013
 * @date last modification: Fri Jan 22 2016
 *
 * @brief  Helper file for Voigt notation
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
template <UInt dim> class VoigtHelper {
public:
  /// transfer the B matrix to a Voigt notation B matrix
  inline static void transferBMatrixToSymVoigtBMatrix(
      const Matrix<Real> & B, Matrix<Real> & Bvoigt, UInt nb_nodes_per_element);

  /// transfer the BNL matrix to a Voigt notation B matrix (See Bathe et al.
  /// IJNME vol 9, 1975)
  inline static void transferBMatrixToBNL(const Matrix<Real> & B,
                                          Matrix<Real> & Bvoigt,
                                          UInt nb_nodes_per_element);

  /// transfer the BL2 matrix to a Voigt notation B matrix (See Bathe et al.
  /// IJNME vol 9, 1975)
  inline static void transferBMatrixToBL2(const Matrix<Real> & B,
                                          const Matrix<Real> & grad_u,
                                          Matrix<Real> & Bvoigt,
                                          UInt nb_nodes_per_element);

public:
  const static UInt size;
  // matrix of vector index I as function of tensor indices i,j
  const static UInt mat[dim][dim];
  // array of matrix indices ij as function of vector index I
  const static UInt vec[dim * dim][2];
  // factors to multiply the strain by for voigt notation
  const static Real factors[dim * (dim - (dim - 1) / 2)];
};

__END_AKANTU__

#include "aka_voigthelper_tmpl.hh"

#endif
