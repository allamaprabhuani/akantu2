/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKA_VOIGTHELPER_HH_
#define AKA_VOIGTHELPER_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim> class AKANTU_EXPORT VoigtHelper {
  static_assert(dim > 0U, "Cannot be < 1D");
  static_assert(dim < 4U, "Cannot be > 3D");

public:
  /* ------------------------------------------------------------------------ */
  template <class M, class V>
  static constexpr inline void matrixToVoigt(M && matrix, V && vector);

  template <class M>
  static constexpr inline decltype(auto) matrixToVoigt(M && matrix);

  template <class M, class V>
  static constexpr inline void matrixToVoigtWithFactors(M && matrix,
                                                        V && vector);

  template <class M>
  static constexpr inline decltype(auto) matrixToVoigtWithFactors(M && matrix);

  template <class M, class V>
  static constexpr inline void voigtToMatrix(V && vector, M && matrix);

  template <class V>
  static constexpr inline decltype(auto) voigtToMatrix(V && vector);
  /* ------------------------------------------------------------------------ */

  /// transfer the B matrix to a Voigt notation B matrix
  template <typename D1, typename D2>
  static constexpr inline void
  transferBMatrixToSymVoigtBMatrix(const Eigen::MatrixBase<D1> & B,
                                   Eigen::MatrixBase<D2> & Bvoigt,
                                   Int nb_nodes_per_element);

  /// transfer the BNL matrix to a Voigt notation B matrix (See Bathe et al.
  /// IJNME vol 9, 1975)
  template <typename D1, typename D2>
  static constexpr inline void
  transferBMatrixToBNL(const Eigen::MatrixBase<D1> & B,
                       Eigen::MatrixBase<D2> & Bvoigt,
                       Int nb_nodes_per_element);

  /// transfer the BL2 matrix to a Voigt notation B matrix (See Bathe et al.
  /// IJNME vol 9, 1975)
  template <typename D1, typename D2, typename D3>
  static constexpr inline void transferBMatrixToBL2(
      const Eigen::MatrixBase<D1> & B, const Eigen::MatrixBase<D2> & grad_u,
      Eigen::MatrixBase<D3> & Bvoigt, Int nb_nodes_per_element);

public:
  static constexpr Int size{(dim * (dim - 1)) / 2 + dim};
  // matrix of vector index I as function of tensor indices i,j
  static const Idx mat[dim][dim];
  // array of matrix indices ij as function of vector index I
  static const Idx vec[dim * dim][2];
  // factors to multiply the strain by for voigt notation
  static const Real factors[size];
};

} // namespace akantu

#include "aka_voigthelper_tmpl.hh"

#endif
