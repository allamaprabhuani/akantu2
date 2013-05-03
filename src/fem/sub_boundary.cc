/**
 * @file   boundary.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 *
 * @date   Wed Mar 06 09:30:00 2013
 *
 * @brief  Stores information relevent to the notion of domain boundary and surfaces.
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
#include <sstream>
#include <algorithm>
#include <iterator>
#include "sub_boundary.hh"
#include "mesh.hh"
#include "aka_csr.hh"
#include "mesh_utils.hh"

__BEGIN_AKANTU__


SubBoundary::SubBoundary(const std::string & _name, const ID & _id, const ID & parent_id, const MemoryID & mem_id)
:Dumpable(parent_id + "_" + _id), id(parent_id + "_" + _id), memory_id(mem_id), name(_name), nodes(0, 1, id + ":nodes"), elements("elements", id, memory_id)
{
}

/* -------------------------------------------------------------------------- */
void SubBoundary::printself(std::ostream & stream) const {
}

/* -------------------------------------------------------------------------- */
void SubBoundary::cleanUpNodeList() {
  UInt size_before = nodes.getSize();
  UInt * newEnd;
  UInt * begin = nodes.storage();
  UInt * end = nodes.storage()+nodes.getSize();
  std::sort(begin, end);
  newEnd = std::unique(begin, end);
  UInt crop_size = end-newEnd;
  UInt size_after = size_before - crop_size;
  nodes.resize(size_after);
}

__END_AKANTU__

