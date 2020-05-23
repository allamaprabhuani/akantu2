/**

 * @file   petsc_wrapper.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Sat Feb 03 2018
 *
 * @brief  Wrapper of PETSc structures
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <petscsys.h>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PETSC_WRAPPER_HH__
#define __AKANTU_PETSC_WRAPPER_HH__

/* -------------------------------------------------------------------------- */
#define PETSc_call(func, ...)                                                  \
  do {                                                                         \
    auto ierr = func(__VA_ARGS__);                                             \
    if (PetscUnlikely(ierr != 0)) {                                            \
      const char * desc;                                                       \
      PetscErrorMessage(ierr, &desc, nullptr);                                 \
      AKANTU_EXCEPTION("Error in PETSc call to \'" << #func                    \
                                                   << "\': " << desc);         \
    }                                                                          \
  } while (false)

namespace akantu {
namespace detail {
  template <typename T> void PETScSetName(T t, const ID & id) {
    PETSc_call(PetscObjectSetName, reinterpret_cast<PetscObject>(t),
               id.c_str());
  }
} // namespace detail
} // namespace akantu

#endif /* __AKANTU_PETSC_WRAPPER_HH__ */
