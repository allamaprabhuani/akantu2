/**
 * @file   aka_plane.hh
 * @author Alejandro Aragon <alejandro.aragon@epfl.ch>
 * @date   Tue Sep 04 15:33:00 2012
 *
 * @brief  class of a plane
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

#ifndef __AKANTU_AKA_PLANE_HH__
#define __AKANTU_AKA_PLANE_HH__

#include "aka_point.hh"


__BEGIN_AKANTU__


struct Plane {
  
  typedef Point<3> point_type;
  typedef point_type::value_type value_type;
  
  Plane(const point_type& a, const point_type& b, const point_type& c) {
    
    n_ = cross(b - a, c - a).normalize();
    d_ = n_*a;
  }
  
  const point_type& normal() const
  { return n_; }
  
  value_type distance() const
  { return d_; }
  
  friend std::ostream& operator<<(std::ostream& os, const Plane& pi) {
    
    os<<"Plane[normal: "<<pi.normal()<<", origin distance: "<<pi.distance()<<"]"<<std::endl;
    return os;
  }
  
private:
  
  point_type n_;     //!< Plane unit normal
  value_type d_;     //!< Distance from the plane to the origin
  
};


__END_AKANTU__

#endif /* __AKANTU_AKA_PLANE_HH__ */
