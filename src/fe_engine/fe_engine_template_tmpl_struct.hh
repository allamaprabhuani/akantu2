/**
 * Copyright (©) 2014-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "shape_structural.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <ElementKind kind, typename = void>
struct AssembleFieldMatrixStructHelper {};

template <ElementKind kind>
struct AssembleFieldMatrixStructHelper<
    kind, typename std::enable_if<kind == _ek_structural>::type> {
  template <template <ElementKind, class> class I,
            template <ElementKind> class S, ElementKind k, class IOF>
  static void call(const FEEngineTemplate<I, S, k, IOF> & fem,
                   const Array<Real> & field_1, Int nb_degree_of_freedom,
                   SparseMatrix & M, Array<Real> * n,
                   ElementTypeMapArray<Real> & rotation_mat, ElementType type,
                   GhostType ghost_type) {
#define ASSEMBLE_MATRIX(type)                                                  \
  fem.template assembleFieldMatrix<type>(field_1, nb_degree_of_freedom, M, n,  \
                                         rotation_mat, ghost_type)

    AKANTU_BOOST_KIND_ELEMENT_SWITCH(ASSEMBLE_MATRIX, _ek_structural);
#undef ASSEMBLE_MATRIX
  }
};

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline void
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::assembleFieldMatrix(
    const Array<Real> &, Int, SparseMatrix &, Array<Real> *,
    ElementTypeMapArray<Real> &, GhostType) const {
  AKANTU_TO_IMPLEMENT();
}

} // namespace akantu
