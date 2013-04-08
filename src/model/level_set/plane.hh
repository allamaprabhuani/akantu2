/**
 * @file   plane.hh
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

#include "aka_common.hh"
#include "geometry.hh"

#ifndef __AKANTU_GEOMETRY_PLANE_HH__
#define	__AKANTU_GEOMETRY_PLANE_HH__

__BEGIN_AKANTU__



class plane : public geometry {
public:
    plane(UInt dim);
    void setNormalPoint(Real * normal, Real * point);
    Real distance(Real * node);
    virtual ~plane();
private:
    Real * Normal;
    Real * Point;
    UInt d;

};

__END_AKANTU__

#endif	/* __AKANTU_GEOMETRY_SPHERE_HH__ */

