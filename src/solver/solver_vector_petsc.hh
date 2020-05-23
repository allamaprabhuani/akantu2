/**
 * @file   solver_vector_petsc.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Jan 01 2019
 *
 * @brief A Documented file.
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
/* -------------------------------------------------------------------------- */
#include "dof_manager_petsc.hh"
#include "vector_petsc.hh"
#include "solver_vector.hh"
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLVER_VECTOR_PETSC_HH__
#define __AKANTU_SOLVER_VECTOR_PETSC_HH__

namespace akantu {

using SolverVectorPETSc = SolverVectorTmpl<VectorPETSc, DOFManagerPETSc>;

} // namespace akantu

#endif /* __AKANTU_SOLVER_VECTOR_PETSC_HH__ */
