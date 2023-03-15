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
#include "aka_voigthelper.hh"
/* -------------------------------------------------------------------------- */

// #ifndef __AKANTU_AKA_VOIGTHELPER_TMPL_HH__
// #define __AKANTU_AKA_VOIGTHELPER_TMPL_HH__

namespace akantu {

template <Int dim> constexpr Int VoigtHelper<dim>::size;

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class M, class V>
constexpr inline void VoigtHelper<dim>::matrixToVoigt(M && matrix, V && vector) {
  for (Int I = 0; I < size; ++I) {
    auto i = vec[I][0];
    auto j = vec[I][1];
    vector(I) = matrix(i, j);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class M>
constexpr inline decltype(auto) VoigtHelper<dim>::matrixToVoigt(M && matrix) {
  Vector<Real, size> vector;
  matrixToVoigt(std::forward<M>(matrix), vector);
  return vector;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class M, class V>
constexpr inline void VoigtHelper<dim>::matrixToVoigtWithFactors(M && matrix,
                                                       V && vector) {
  for (Int I = 0; I < size; ++I) {
    auto i = vec[I][0];
    auto j = vec[I][1];
    vector(I) = factors[I] * matrix(i, j);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class M>
constexpr inline decltype(auto) VoigtHelper<dim>::matrixToVoigtWithFactors(M && matrix) {
  Vector<Real, size> vector;
  matrixToVoigtWithFactors(std::forward<M>(matrix), vector);
  return vector;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class M, class V>
constexpr inline void VoigtHelper<dim>::voigtToMatrix(V && vector, M && matrix) {
  for (Int I = 0; I < size; ++I) {
    auto i = vec[I][0];
    auto j = vec[I][1];
    matrix(i, j) = matrix(j, i) = vector(I);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class V>
constexpr inline decltype(auto) VoigtHelper<dim>::voigtToMatrix(V && vector) {
  Matrix<Real, dim, dim> matrix;
  voigtToMatrix(std::forward<V>(vector), matrix);
  return matrix;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename D1, typename D2>
constexpr inline void VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(
    const Eigen::MatrixBase<D1> & B, Eigen::MatrixBase<D2> & Bvoigt,
    Int nb_nodes_per_element) {
  Bvoigt.zero();

  for (Int i = 0; i < dim; ++i) {
    for (Int n = 0; n < nb_nodes_per_element; ++n) {
      Bvoigt(i, i + n * dim) = B(i, n);
    }
  }

  if (dim == 2) {
    /// in 2D, fill the @f$ [\frac{\partial N_i}{\partial x}, \frac{\partial
    /// N_i}{\partial y}]@f$ row
    for (Int n = 0; n < nb_nodes_per_element; ++n) {
      Bvoigt(2, 1 + n * 2) = B(0, n);
      Bvoigt(2, 0 + n * 2) = B(1, n);
    }
  }

  if (dim == 3) {
    for (Int n = 0; n < nb_nodes_per_element; ++n) {
      auto dndx = B(0, n);
      auto dndy = B(1, n);
      auto dndz = B(2, n);

      /// in 3D, fill the @f$ [0, \frac{\partial N_i}{\partial y},
      /// \frac{N_i}{\partial z}]@f$ row
      Bvoigt(3, 1 + n * 3) = dndz;
      Bvoigt(3, 2 + n * 3) = dndy;

      /// in 3D, fill the @f$ [\frac{\partial N_i}{\partial x}, 0,
      /// \frac{N_i}{\partial z}]@f$ row
      Bvoigt(4, 0 + n * 3) = dndz;
      Bvoigt(4, 2 + n * 3) = dndx;

      /// in 3D, fill the @f$ [\frac{\partial N_i}{\partial x},
      /// \frac{N_i}{\partial y}, 0]@f$ row
      Bvoigt(5, 0 + n * 3) = dndy;
      Bvoigt(5, 1 + n * 3) = dndx;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename D1, typename D2>
constexpr inline void
VoigtHelper<dim>::transferBMatrixToBNL(const Eigen::MatrixBase<D1> & B,
                                       Eigen::MatrixBase<D2> & Bvoigt,
                                       Int nb_nodes_per_element) {
  Bvoigt.zero();

  // see Finite element formulations for large deformation dynamic analysis,
  // Bathe et al. IJNME vol 9, 1975, page 364 B_{NL}
  for (Int i = 0; i < dim; ++i) {
    for (Int m = 0; m < nb_nodes_per_element; ++m) {
      for (Int n = 0; n < dim; ++n) {
        // std::cout << B(n, m) << std::endl;
        Bvoigt(i * dim + n, m * dim + i) = B(n, m);
      }
    }
  }
  // TODO: Verify the 2D and 1D case
}

/* -------------------------------------------------------------------------- */
template <>
template <typename D1, typename D2, typename D3>
constexpr inline void VoigtHelper<1>::transferBMatrixToBL2(
    const Eigen::MatrixBase<D1> & B, const Eigen::MatrixBase<D2> & grad_u,
    Eigen::MatrixBase<D3> & Bvoigt, Int nb_nodes_per_element) {
  Bvoigt.zero();
  for (Int j = 0; j < nb_nodes_per_element; ++j) {
    Bvoigt(0, j) = grad_u(0, 0) * B(0, j);
  }
}

/* -------------------------------------------------------------------------- */
template <>
template <typename D1, typename D2, typename D3>
constexpr inline void VoigtHelper<3>::transferBMatrixToBL2(
    const Eigen::MatrixBase<D1> & dNdX, const Eigen::MatrixBase<D2> & grad_u,
    Eigen::MatrixBase<D3> & Bvoigt, Int nb_nodes_per_element) {
  Bvoigt.zero();

  for (Int I = 0; I < 3; ++I) {
    for (Int a = 0; a < nb_nodes_per_element; ++a) {
      for (Int i = 0; i < 3; ++i) {
        Bvoigt(I, a * 3 + i) = grad_u(i, I) * dNdX(I, a);
      }
    }
  }

  for (Int Iv = 3; Iv < 6; ++Iv) {
    for (Int a = 0; a < nb_nodes_per_element; ++a) {
      for (Int k = 0; k < 3; ++k) {
        auto aux = Iv - 3;
        for (Int m = 0; m < 3; ++m) {
          if (m != aux) {
            auto index1 = m;
            auto index2 = 3 - m - aux;
            Bvoigt(Iv, a * 3 + k) += grad_u(k, index1) * dNdX(index2, a);
          }
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <>
template <typename D1, typename D2, typename D3>
constexpr inline void VoigtHelper<2>::transferBMatrixToBL2(
    const Eigen::MatrixBase<D1> & B, const Eigen::MatrixBase<D2> & grad_u,
    Eigen::MatrixBase<D3> & Bvoigt, Int nb_nodes_per_element) {

  Bvoigt.zero();

  for (Int i = 0; i < 2; ++i) {
    for (Int j = 0; j < nb_nodes_per_element; ++j) {
      for (Int k = 0; k < 2; ++k) {
        Bvoigt(i, j * 2 + k) = grad_u(k, i) * B(i, j);
      }
    }
  }

  for (Int j = 0; j < nb_nodes_per_element; ++j) {
    for (Int k = 0; k < 2; ++k) {
      for (Int m = 0; m < 2; ++m) {
        auto index1 = m;
        auto index2 = (2 - 1) - m;
        Bvoigt(2, j * 2 + k) += grad_u(k, index1) * B(index2, j);
      }
    }
  }
}

} // namespace akantu

//#endif /* __AKANTU_AKA_VOIGTHELPER_TMPL_HH__ */
