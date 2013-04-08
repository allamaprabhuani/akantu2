/**
 * @file   sphere.cc
 *
 * @author Danie Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Fri Dec 14 11:45:43 2012
 *
 * @brief  Geometry of a sphere
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
#include "sphere.hh"

__BEGIN_AKANTU__

sphere::sphere(UInt dim) {
    d = dim;
    if (d)
        Center = new Real [d];
}

void sphere::setCenterRadius(Real * center, Real r) {
    for (UInt i = 0; i < d; i++)
        Center[i] = center[i];
    R = r;
    
    //AKANTU_DEBUG_ASSERT(R > Math::getTolerance(), "The Radius is wrong");
}

Real sphere::distance(Real * node) {
    Real dist = 0.0;
    for (UInt i = 0; i < d; i++)
        dist += (node[i] - Center[i])*(node[i] - Center[i]);
    dist = sqrt(dist)-R;
    return dist;
}

sphere::~sphere() {
    delete [] Center;
}

__END_AKANTU__