/**
 * @file   container.cc
 *
 * @author Danie Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Fri Dec 14 11:45:43 2012
 *
 * @brief  Geometry of a container
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

#include "aka_math.hh"
#include "container.hh"
#include "geometry.hh"
#include "mesh.hh"

__BEGIN_AKANTU__

container::container() {
    ;
}

Real container::distance(Real * node) {
    UInt n_geos = Geometries.getSize();
    if (n_geos) {
        geometry * geo = *(Geometries.storage());
        Real dist_min = geo->distance(node);
        for (UInt i = 1; i < n_geos; i++) {
            geo=*(Geometries.storage()+i);
            Real dist = geo->distance(node);
            if (dist < dist_min)
                dist_min = dist;
        }
        return dist_min;
    } else
        AKANTU_DEBUG_ERROR("Container empty!");
}

container::~container() {
    ;
}

__END_AKANTU__