/**
 * Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "contact_element.hh"
#include "contact_mechanics_model.hh"
#include "fe_engine.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_RESOLUTION_UTILS_HH__
#define __AKANTU_RESOLUTION_UTILS_HH__

/* -------------------------------------------------------------------------- */

namespace akantu {

class ResolutionUtils {
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// computes the shape function matric for the contact element (@f$A
  /// @f$) where row is equal to spatial dimension and cols is equal
  /// to spatial dimension times number of nodes in contact element
  template <class Derived>
  static Matrix<Real>
  computeShapeFunctionMatrix(const ContactElement & element,
                             const Eigen::MatrixBase<Derived> & projection) {
    const auto type = element.master.type;
    const auto surface_dimension = Mesh::getSpatialDimension(type);
    const auto spatial_dimension = surface_dimension + 1;
    const auto nb_nodes_per_contact = element.getNbNodes();

    Matrix<Real> shape_matrix(spatial_dimension,
                              spatial_dimension * nb_nodes_per_contact);
    shape_matrix.zero();

    auto && shapes = ElementClassHelper<_ek_regular>::getN(projection, type);

    for (auto i : arange(nb_nodes_per_contact)) {
      for (auto j : arange(spatial_dimension)) {
        if (i == 0) {
          shape_matrix(j, i * spatial_dimension + j) = 1;
          continue;
        }
        shape_matrix(j, i * spatial_dimension + j) = -shapes[i - 1];
      }
    }
    return shape_matrix;
  }

  template <class Derived>
  static Tensor3<Real> computeDerivativeShapeFunctionMatrix(
      const ContactElement & element,
      const Eigen::MatrixBase<Derived> & projection) {

    const auto type = element.master.type;

    // computing shape derivatives
    auto && shape_derivatives =
        ElementClassHelper<_ek_regular>::getDNDS(projection, type);

    const auto surface_dimension = shape_derivatives.rows();
    const auto spatial_dimension = shape_derivatives.rows() + 1;
    const auto nb_nodes_per_contact = shape_derivatives.cols() + 1;

    Tensor3<Real> derivative_shape_matrix(
        spatial_dimension, spatial_dimension * nb_nodes_per_contact,
        surface_dimension);
    derivative_shape_matrix.zero();

    for (auto && [dnds, Aj] :
         zip(shape_derivatives.transpose(), derivative_shape_matrix)) {
      for (auto i : arange(nb_nodes_per_contact)) {
        for (auto j : arange(spatial_dimension)) {
          if (i == 0) {
            Aj(j, i * spatial_dimension + j) = 0;
            continue;
          }
          Aj(j, i * spatial_dimension + j) = dnds(i - 1);
        }
      }
    }
    return derivative_shape_matrix;
  }
};

} // namespace akantu

#endif /* __AKANTU_RESOLUTION_UTILS_HH__ */
