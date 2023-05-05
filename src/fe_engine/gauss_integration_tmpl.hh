/**
 * Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_math.hh"

#ifndef AKANTU_GAUSS_INTEGRATION_TMPL_HH_
#define AKANTU_GAUSS_INTEGRATION_TMPL_HH_

#include "element_class.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
/* GaussIntegrationElement                                                    */
/* -------------------------------------------------------------------------- */
namespace _aka_gauss_helpers {
  template <GaussIntegrationType type, Int n>
  struct GaussIntegrationNbPoints {};

  template <Int n> struct GaussIntegrationNbPoints<_git_not_defined, n> {
    static constexpr Int nb_points = 0;
  };

  template <Int n> struct GaussIntegrationNbPoints<_git_point, n> {
    static constexpr Int nb_points = 1;
  };

  template <Int n> struct GaussIntegrationNbPoints<_git_segment, n> {
    static constexpr Int nb_points = (n + 1) / 2 + (bool((n + 1) % 2) ? 1 : 0);
  };

#define DECLARE_GAUSS_NB_POINTS(type, order, points)                           \
  template <> struct GaussIntegrationNbPoints<type, order> {                   \
    static constexpr Int nb_points = points;                                   \
  }

#define DECLARE_GAUSS_NB_POINTS_PENT(type, order, xo, yo)                      \
  template <> struct GaussIntegrationNbPoints<type, order> {                   \
    static constexpr Int x_order = xo;                                         \
    static constexpr Int yz_order = yo;                                        \
    static constexpr Int nb_points = 1;                                        \
  }

  DECLARE_GAUSS_NB_POINTS(_git_triangle, 1, 1);
  DECLARE_GAUSS_NB_POINTS(_git_triangle, 2, 3);
  DECLARE_GAUSS_NB_POINTS(_git_triangle, 3, 4);
  DECLARE_GAUSS_NB_POINTS(_git_triangle, 4, 6);
  DECLARE_GAUSS_NB_POINTS(_git_triangle, 5, 7);
  DECLARE_GAUSS_NB_POINTS(_git_tetrahedron, 1, 1);
  DECLARE_GAUSS_NB_POINTS(_git_tetrahedron, 2, 4);
  DECLARE_GAUSS_NB_POINTS(_git_tetrahedron, 3, 5);
  DECLARE_GAUSS_NB_POINTS(_git_tetrahedron, 4, 15);
  DECLARE_GAUSS_NB_POINTS(_git_tetrahedron, 5, 15);
  DECLARE_GAUSS_NB_POINTS_PENT(_git_pentahedron, 1, 3,
                               2); // order 3 in x, order 2 in y and z
  DECLARE_GAUSS_NB_POINTS_PENT(_git_pentahedron, 2, 3,
                               2); // order 3 in x, order 2 in y and z
  DECLARE_GAUSS_NB_POINTS_PENT(_git_pentahedron, 3, 3,
                               3); // order 3 in x, order 3 in y and z
  DECLARE_GAUSS_NB_POINTS_PENT(_git_pentahedron, 4, 5,
                               5); // order 5 in x, order 5 in y and z
  DECLARE_GAUSS_NB_POINTS_PENT(_git_pentahedron, 5, 5,
                               5); // order 5 in x, order 5 in y and z

  template <GaussIntegrationType type, Int n, Int on = n,
            bool end_recurse = false>
  struct GaussIntegrationNbPointsHelper {
    static constexpr Int pnp = GaussIntegrationNbPoints<type, n>::nb_points;
    static constexpr Int order = n;
    static constexpr Int nb_points = pnp;
  };

  template <GaussIntegrationType type, Int n, Int on>
  struct GaussIntegrationNbPointsHelper<type, n, on, true> {
    static constexpr Int nb_points = 0;
  };

  /* ------------------------------------------------------------------------ */
  /* Generic helper                                                           */
  /* ------------------------------------------------------------------------ */
  template <GaussIntegrationType type, Int dimension, Int n>
  struct GaussIntegrationTypeDataHelper {
    using git_np = GaussIntegrationNbPoints<type, n>;
    using git_data = GaussIntegrationTypeData<type, git_np::nb_points>;

    static constexpr Int getNbQuadraturePoints() { return git_np::nb_points; }

    static constexpr auto getQuadraturePoints()
        -> Matrix<Real, std::max(1, dimension),
                  GaussIntegrationTypeDataHelper::getNbQuadraturePoints()> {
      constexpr auto nb_points = git_np::nb_points;
      using Matrix = Eigen::Matrix<Real, std::max(1, dimension), nb_points>;
      Matrix quads = Eigen::Map<const Matrix>(git_data::quad_positions);
      return quads;
    }

    static constexpr auto getWeights()
        -> Vector<Real,
                  GaussIntegrationTypeDataHelper::getNbQuadraturePoints()> {
      constexpr auto nb_points = git_np::nb_points;
      using Vector = Eigen::Matrix<Real, nb_points, 1>;
      Vector weights = Eigen::Map<const Vector>(git_data::quad_weights);
      return weights;
    }
  };

#if !defined(DOXYGEN)
  /* ------------------------------------------------------------------------ */
  /* helper for _segment _quadrangle _hexahedron                              */
  /* ------------------------------------------------------------------------ */
  template <Int dimension, Int dp>
  struct GaussIntegrationTypeDataHelper<_git_segment, dimension, dp> {
    using git_np = GaussIntegrationNbPoints<_git_segment, dp>;
    using git_data = GaussIntegrationTypeData<_git_segment, git_np::nb_points>;

    static constexpr Int getNbQuadraturePoints() {
      return Math::pow(git_np::nb_points, dimension);
    }

    static constexpr auto getQuadraturePoints()
        -> Matrix<Real, dimension,
                  GaussIntegrationTypeDataHelper::getNbQuadraturePoints()> {
      const auto tot_nquad = getNbQuadraturePoints();
      const auto nquad = git_np::nb_points;

      Eigen::Matrix<Real, dimension, tot_nquad> quads;
      Eigen::Map<const Eigen::Matrix<Real, nquad, 1>> pos(
          git_data::quad_positions);

      Int offset = 1;
      for (Int d = 0; d < dimension; ++d) {
        for (Int n = 0, q = 0; n < tot_nquad; ++n, q += offset) {
          Int rq = q % tot_nquad + q / tot_nquad;
          quads(d, rq) = pos(n % nquad);
        }
        offset *= nquad;
      }
      return quads;
    }

    static constexpr Vector<
        Real, GaussIntegrationTypeDataHelper::getNbQuadraturePoints()>
    getWeights() {
      const auto tot_nquad = getNbQuadraturePoints();
      const auto nquad = git_np::nb_points;

      Eigen::Matrix<Real, tot_nquad, 1> quads_weights;
      quads_weights.fill(1);
      Eigen::Map<const Eigen::Matrix<Real, nquad, 1>> weights(
          git_data::quad_weights);

      Int offset = 1;
      for (Int d = 0; d < dimension; ++d) {
        for (Int n = 0, q = 0; n < tot_nquad; ++n, q += offset) {
          Int rq = q % tot_nquad + q / tot_nquad;
          quads_weights(rq) *= weights(n % nquad);
        }
        offset *= nquad;
      }
      return quads_weights;
    }
  };

  /* ------------------------------------------------------------------------ */
  /* helper for _pentahedron                                                  */
  /* ------------------------------------------------------------------------ */
  template <Int dimension, Int dp>
  struct GaussIntegrationTypeDataHelper<_git_pentahedron, dimension, dp> {
    using git_info = GaussIntegrationNbPoints<_git_pentahedron, dp>;
    using git_np_seg =
        GaussIntegrationNbPoints<_git_segment, git_info::x_order>;
    using git_np_tri =
        GaussIntegrationNbPoints<_git_triangle, git_info::yz_order>;
    using git_data_seg =
        GaussIntegrationTypeData<_git_segment, git_np_seg::nb_points>;
    using git_data_tri =
        GaussIntegrationTypeData<_git_triangle, git_np_tri::nb_points>;

    static constexpr Int getNbQuadraturePoints() {
      return git_np_seg::nb_points * git_np_tri::nb_points;
    }

    static constexpr auto getQuadraturePoints()
        -> Matrix<Real, dimension,
                  GaussIntegrationTypeDataHelper::getNbQuadraturePoints()> {
      constexpr auto tot_nquad = getNbQuadraturePoints();
      constexpr auto nquad_seg = git_np_seg::nb_points;
      constexpr auto nquad_tri = git_np_tri::nb_points;

      Matrix<Real, dimension, tot_nquad> quads;
      Eigen::Map<Vector<Real, nquad_seg>> pos_seg(git_data_seg::quad_positions);
      Eigen::Map<Matrix<Real, 2, nquad_tri>> pos_tri(
          git_data_tri::quad_positions);

      for (Int ns = 0, q = 0; ns < nquad_seg; ++ns) {
        for (Int nt = 0; nt < nquad_tri; ++nt, ++q) {
          auto && quad = quads(q);
          quad(_x) = pos_seg(ns);
          quad(_y) = pos_tri(_x, nt);
          quad(_z) = pos_tri(_y, nt);
        }
      }
      return quads;
    }

    static constexpr auto getWeights()
        -> Vector<Real,
                  GaussIntegrationTypeDataHelper::getNbQuadraturePoints()> {
      constexpr auto tot_nquad = getNbQuadraturePoints();
      constexpr auto nquad_seg = git_np_seg::nb_points;
      constexpr auto nquad_tri = git_np_tri::nb_points;

      Vector<Real, tot_nquad> quads_weights;
      Eigen::Map<const Vector<Real, nquad_seg>> weight_seg(
          git_data_seg::quad_weights);
      Eigen::Map<const Vector<Real, nquad_tri>> weight_tri(
          git_data_tri::quad_weights);

      for (Int ns = 0, q = 0; ns < nquad_seg; ++ns) {
        for (Int nt = 0; nt < nquad_tri; ++nt, ++q) {
          quads_weights(q) = weight_seg(ns) * weight_tri(nt);
        }
      }
      return quads_weights;
    }
  };
#endif
} // namespace _aka_gauss_helpers

/* -------------------------------------------------------------------------- */
template <ElementType element_type, Int n>
constexpr Int
GaussIntegrationElement<element_type, n>::getNbQuadraturePoints() {
  using data_helper = _aka_gauss_helpers::GaussIntegrationTypeDataHelper<
      ElementClassProperty<element_type>::gauss_integration_type,
      interpolation_property::natural_space_dimension, n>;
  return data_helper::getNbQuadraturePoints();
}

/* -------------------------------------------------------------------------- */
template <ElementType element_type, Int n>
constexpr auto GaussIntegrationElement<element_type, n>::getWeights() {
  using data_helper = _aka_gauss_helpers::GaussIntegrationTypeDataHelper<
      ElementClassProperty<element_type>::gauss_integration_type,
      interpolation_property::natural_space_dimension, n>;
  return data_helper::getWeights();
}

template <ElementType type, Int n>
constexpr auto GaussIntegrationElement<type, n>::getQuadraturePoints() {
  using data_helper = _aka_gauss_helpers::GaussIntegrationTypeDataHelper<
      ElementClassProperty<type>::gauss_integration_type,
      interpolation_property::natural_space_dimension, n>;
  return data_helper::getQuadraturePoints();
}

} // namespace akantu

#endif /* AKANTU_GAUSS_INTEGRATION_TMPL_HH_ */
