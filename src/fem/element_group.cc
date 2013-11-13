/**
 * @file   element_group.cc
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
#include "mesh.hh"
#include "aka_csr.hh"
#include "mesh_utils.hh"

#include "element_group.hh"
#if defined(AKANTU_USE_IOHELPER)
#  include "dumper_paraview.hh"
#endif

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
ElementGroup::ElementGroup(const std::string & boundary_name,
                           const Mesh & mesh,
                           NodeGroup & node_group,
                           UInt dimension,
                           const std::string & id,
                           const MemoryID & mem_id) :
  Memory(id, mem_id),
  mesh(mesh),
  name(boundary_name),
  elements("elements", id, mem_id),
  node_group(node_group),
  dimension(dimension) {
  AKANTU_DEBUG_IN();

#if defined(AKANTU_USE_IOHELPER)
  this->registerDumper<DumperParaview>("paraview_"  + boundary_name, boundary_name, true);
  this->addDumpFilteredMesh(mesh, elements, node_group.getNodes(), dimension);
#endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ElementGroup::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "ElementGroup [" << std::endl;
  stream << space << " + name: " << name << std::endl;
  elements.printself(stream, indent + 1);
  node_group.printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}

__END_AKANTU__

