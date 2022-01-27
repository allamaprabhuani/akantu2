/**
 * @file   aka_math.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Peter Spijker <peter.spijker@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Tue Sep 29 2020
 *
 * @brief  Implementation of the math toolbox
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_math.hh"
#include "aka_array.hh"
#include "aka_iterators.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
namespace Math {
  /* ---------------------------------------------------------------------- */
  void compute_tangents(const Array<Real> & normals, Array<Real> & tangents) {
    if (normals.empty()) {
      return;
    }

    auto spatial_dimension = normals.getNbComponent();
    auto tangent_components = spatial_dimension * (spatial_dimension - 1);

    if (tangent_components == 0) {
      return;
    }

    AKANTU_DEBUG_ASSERT(
        tangent_components == tangents.getNbComponent(),
        "Cannot compute the tangents, the storage array for tangents"
            << " does not have the good amount of components.");

    auto nb_normals = normals.size();
    tangents.resize(nb_normals);
    tangents.zero();

    /// compute first tangent
    for (auto && data : zip(make_view(normals, spatial_dimension),
                            make_view(tangents, tangent_components))) {
      const auto & normal_ = std::get<0>(data);
      auto & tangent = std::get<1>(data);

      if (are_float_equal(normal_.norm(), 0.)) {
        tangent(0) = 1.;
      } else {
        tangent = normal(normal_);
      }
    }

    /// compute second tangent (3D case)
    if (spatial_dimension == 3) {
      for (auto && data : zip(make_view(normals, spatial_dimension),
                              make_view(tangents, tangent_components))) {
        const auto & normal_ = std::get<0>(data);
        auto & tangent = std::get<1>(data);

        tangent.col(1) = normal(normal_, tangent.col(0));
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  Real reduce(Array<Real> & array) {
    UInt nb_values = array.size();
    if (nb_values == 0) {
      return 0.;
    }

    Int nb_values_to_sum = nb_values >> 1;

    std::sort(array.begin(), array.end());

    // as long as the half is not empty
    while (nb_values_to_sum != 0) {
      Int remaining = (nb_values - 2 * nb_values_to_sum);
      if (remaining != 0U) {
        array(nb_values - 2) += array(nb_values - 1);
      }

      // sum to consecutive values and store the sum in the first half
      for (Int i = 0; i < nb_values_to_sum; ++i) {
        array(i) = array(2 * i) + array(2 * i + 1);
      }

      nb_values = nb_values_to_sum;
      nb_values_to_sum >>= 1;
    }

    return array(0);
  }
} // namespace Math
} // namespace akantu
