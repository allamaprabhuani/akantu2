/**
 * @file   plane.cc
 *
 * @author Danie Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Fri Dec 14 11:45:43 2012
 *
 * @brief  Geometry of a plane
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

#include "plane.hh"
#include "aka_math.hh"

__BEGIN_AKANTU__

plane::plane(UInt dim) {
    d = dim;
    if (d) {
        Normal = new Real [d];
        Point = new Real [d];
    }
}

void plane::setNormalPoint(Real * normal, Real * point) {
    Real norm = 0.0;
    for (UInt i = 0; i < d; i++) {
        Normal[i] = normal[i];
        norm += Normal[i] * Normal[i];
        Point[i] = point[i];
    }
    
    norm = sqrt(norm);
    
    AKANTU_DEBUG_ASSERT(norm > Math::getTolerance(), "The normal is wrong");

    for (UInt i = 0; i < d; i++)
        Normal[i] /= norm;

}

Real plane::distance(Real * node) {
    Real dist = 0.0;
    for (UInt i = 0; i < d; i++)
        dist += (node[i] - Point[i]) * Normal[i];
    return dist;
}

plane::~plane() {
    delete [] Normal;
    delete [] Point;
}

__END_AKANTU__