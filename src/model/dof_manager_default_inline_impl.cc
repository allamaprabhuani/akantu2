/**
 * @file   dof_manager_default_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 12 10:57:47 2015
 *
 * @brief  Implementation of the DOFManagerDefault inline functions
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
#include "dof_manager_default.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_DOF_MANAGER_DEFAULT_INLINE_IMPL_CC__
#define __AKANTU_DOF_MANAGER_DEFAULT_INLINE_IMPL_CC__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline bool DOFManagerDefault::isLocalOrMasterDOF(UInt dof_num) {
  Int dof_type = this->dofs_type(dof_num);
  return (dof_type == -2) || (dof_type == -1);
}

/* -------------------------------------------------------------------------- */
inline const Array<UInt> &
DOFManagerDefault::getLocalEquationNumbers(const ID & dof_id) const {
  return this->getDOFData(dof_id).local_equation_number;
}

inline const Array<UInt> &
DOFManagerDefault::getDOFsAssociatedNodes(const ID & dof_id) const {
  const DOFDataDefault & dof_data = this->getDOFDataTyped<DOFDataDefault>(dof_id);
  return dof_data.associated_nodes;
}
/* -------------------------------------------------------------------------- */
inline void DOFManagerDefault::extractElementEquationNumber(
    const Array<UInt> & equation_numbers, const Vector<UInt> & connectivity,
    UInt nb_degree_of_freedom, Vector<UInt> & element_equation_number) {
  for (UInt i = 0, ld = 0; i < connectivity.size(); ++i) {
    UInt n = connectivity(i);
    for (UInt d = 0; d < nb_degree_of_freedom; ++d, ++ld) {
      element_equation_number(ld) =
          equation_numbers(n * nb_degree_of_freedom + d);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class S>
inline void DOFManagerDefault::localToGlobalEquationNumber(S & inout) {
  for (UInt i = 0; i < inout.size(); ++i) {
    inout(i) = this->global_equation_number(inout(i));
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void DOFManagerDefault::localToGlobalEquationNumber<UInt>(UInt & inout) {
  inout = this->global_equation_number(inout);
}

/* -------------------------------------------------------------------------- */
inline UInt DOFManagerDefault::globalToLocalEquationNumber(UInt global) const {
  equation_numbers_map::const_iterator it =
      this->global_to_local_mapping.find(global);
  AKANTU_DEBUG_ASSERT(it != this->global_to_local_mapping.end(),
                      "This global equation number "
                          << global << " does not exists in " << this->id);

  return it->second;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__

#endif /* __AKANTU_DOF_MANAGER_DEFAULT_INLINE_IMPL_CC__ */
