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
#include "ntn_base_contact.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
inline Int NTNBaseContact::getNbData(const Array<Element> & elements,
                                     const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  Int size = 0;
  Int spatial_dimension = this->model.getSpatialDimension();

  Int nb_nodes = 0;
  for (const auto & el : elements) {
    nb_nodes += Mesh::getNbNodesPerElement(el.type);
  }

  switch (tag) {
  case SynchronizationTag::_cf_nodal: {
    size += nb_nodes * spatial_dimension * sizeof(Real) *
            3; // disp, vel and cur_pos
    break;
  }
  case SynchronizationTag::_cf_incr: {
    size += nb_nodes * spatial_dimension * sizeof(Real) * 1;
    break;
  }
  default: {
  }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void NTNBaseContact::packData(CommunicationBuffer & buffer,
                                     const Array<Element> & elements,
                                     const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  switch (tag) {

  case SynchronizationTag::_cf_nodal: {
    DataAccessor::packNodalDataHelper(this->model.getDisplacement(), buffer,
                                      elements, this->model.getMesh());
    DataAccessor::packNodalDataHelper(this->model.getCurrentPosition(), buffer,
                                      elements, this->model.getMesh());
    DataAccessor::packNodalDataHelper(this->model.getVelocity(), buffer,
                                      elements, this->model.getMesh());
    break;
  }
  case SynchronizationTag::_cf_incr: {
    DataAccessor::packNodalDataHelper(this->model.getIncrement(), buffer,
                                      elements, this->model.getMesh());
    break;
  }
  default: {
  }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void NTNBaseContact::unpackData(CommunicationBuffer & buffer,
                                       const Array<Element> & elements,
                                       const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  switch (tag) {

  case SynchronizationTag::_cf_nodal: {
    DataAccessor::unpackNodalDataHelper(this->model.getDisplacement(), buffer,
                                        elements, this->model.getMesh());
    DataAccessor::unpackNodalDataHelper(
        const_cast<Array<Real> &>(this->model.getCurrentPosition()), buffer,
        elements, this->model.getMesh());
    DataAccessor::unpackNodalDataHelper(this->model.getVelocity(), buffer,
                                        elements, this->model.getMesh());
    break;
  }
  case SynchronizationTag::_cf_incr: {
    DataAccessor::unpackNodalDataHelper(this->model.getIncrement(), buffer,
                                        elements, this->model.getMesh());
    break;
  }
  default: {
  }
  }

  AKANTU_DEBUG_OUT();
}
} // namespace akantu
