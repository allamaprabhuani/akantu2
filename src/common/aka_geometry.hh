/**
 * @file   aka_geometry.hh
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @date   Tue Oct 23 10:35:00 2012
 *
 * @brief  geometric operations
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


#ifndef __AKANTU_GEOMETRY_HH__
#define __AKANTU_GEOMETRY_HH__

#include <iostream>

#include "aka_point.hh"
#include "aka_plane.hh"


__BEGIN_AKANTU__


// predicates


Real left_turn(const Point<2>& p, const Point<2>& q, const Point<2>& r);


//! Tests if point P lies inside a triangle
/*! The triangle is defined by points \c a, \c b and \c c.
 * Note that the triangle points are passed by value because
 * a copy of them is needed.
 */
template <typename T>
bool is_point_in_triangle(const Point<3,T>& p,
                          Point<3,T> a,
                          Point<3,T> b,
                          Point<3,T> c) {
  
  typedef Point<3,T> point_type;
  
  // translate point and triangle so that point lies at origin
  a -= p; b -= p; c -= p;
  
  // compute normal vectors for triangles pab and pbc
  point_type u = cross(b, c);
  point_type v = cross(c, a);
  
  // make sure they are both pointing in the same direction
  if ((u*v) < 0.)
    return false;
  
  // compute normal vector for triangle pca
  point_type w = cross(a, b);
  
  // make sure it points in the same direction as the first two
  if (u*w < 0.)
    return false;
  
  // otherwise P must be in (or on) the triangle
  return true;
}


//! Tests if point P lies inside a triangle
/*! The triangle is defined by points \c a, \c b and \c c.
 * This is a cheaper version of the function \c is_point_in_triangle.
 * Note that the triangle points are passed by value because
 * a copy of them is needed.
 */
template <typename T>
bool is_point_in_triangle2(const Point<3,T>& p,
                           Point<3,T> a,
                           Point<3,T> b,
                           Point<3,T> c) {
  
  // translate point and triangle so that point lies at origin
  a -= p; b -= p; c -= p;
  T ab = a*b;
  T ac = a*c;
  T bc = b*c;
  T cc = c*c;
  
  // make sure plane normals for pab and pca point in the same direction
  if (bc*ac - cc*ab < 0.)
    return false;
  
  // make sure plane normals for pab and pbc point in the same direction
  T bb = b*b;
  if (ab*bc - ac*bb < 0.)
    return false;
  
  // otherwiese p must lie in (or on) the triangle
  return true;
}

// closest point computations


//! Computes the closest point laying on a segment to a point
/*! Given segment \c ab and point \c c, computes closest point \c d on ab.
 * Also returns \c t for the position of the point: a + t*(b - a)
 */
template <typename T>
Point<3,T> closest_point_to_segment(const Point<3,T>& c,
                                    const Point<3,T>& a,
                                    const Point<3,T>& b) {
  
  Point<3,T> ab = b - a;
  
  
  // project c onto ab, computing parameterized position d(t) = a + t*(b – a)
  
  T t = (c - a)*ab / (ab*ab);
  
  // if outside segment, clamp t (and therefore d) to the closest endpoint
  if (t < 0.)
    t = 0.;
  if (t > 1.)
    t = 1.;
  
  // compute projected position from the clamped t
  return a + t * ab;
}



//! Compute the closest point to a triangle
/*! This function uses the concept of Voronoi regions to determine
 * the closest point \c p to a triangle defined by points \c a, \c b
 * \c c.
 */
template <typename T>
Point<3,T> closest_point_to_triangle(const Point<3,T>& p,
                                     const Point<3,T>& a,
                                     const Point<3,T>& b,
                                     const Point<3,T>& c) {
  
  typedef Point<3,T> point_type;
  
  // check if P in vertex region outside A
  point_type ab = b - a;
  point_type ac = c - a;
  point_type ap = p - a;
  
  // compute scalar products
  T d1 = ab * ap;
  T d2 = ac * ap;
  
  if (d1 <= 0. && d2 <= 0.)
    return a; // barycentric coordinates (1,0,0)
  
  // check if P in vertex region outside B
  point_type bp = p - b;
  
  T d3 = ab * bp;
  T d4 = ac * bp;
  
  if (d3 >= 0.0f && d4 <= d3)
    return b; // barycentric coordinates (0,1,0)
  
  // check if P in edge region of AB, if so return projection of P onto AB
  T vc = d1*d4 - d3*d2;
  if (vc <= 0. && d1 >= 0. && d3 <= 0.) {
    T v = d1 / (d1 - d3);
    return a + v * ab; // barycentric coordinates (1-v,v,0)
  }
  
  // check if P in vertex region outside C
  point_type cp = p - c;
  T d5 = ab * cp;
  T d6 = ac * cp;
  if (d6 >= 0.0f && d5 <= d6)
    return c; // barycentric coordinates (0,0,1)
  
  // check if P in edge region of AC, if so return projection of P onto AC
  T vb = d5*d2 - d1*d6;
  if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
    T w = d2 / (d2 - d6);
    return a + w * ac; // barycentric coordinates (1-w,0,w)
  }
  
  // Check if P in edge region of BC, if so return projection of P onto BC
  T va = d3*d6 - d5*d4;
  if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
    T w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return b + w * (c - b); // barycentric coordinates (0,1-w,w)
  }
  
  // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
  T denom = 1.0f / (va + vb + vc);
  T v = vb * denom;
  T w = vc * denom;
  
  return a + ab*v + ac*w; // = u*a + v*b + w*c, u = va*denom = 1.0f - v - w
}


template <typename T>
Point<3,T> closest_ponint_to_plane(const Point<3,T>& q, const Plane& p) {
  
  typedef Point<3,T> point_type;
  
  const point_type& n = p.normal();
  
  T t = (n*q - p.distance()) / (n*n);
  return q - t * n;
}

//! Compute the closest point to a triangle
/*! This function uses the concept of Voronoi regions to determine
 * the closest point \c p to a triangle defined by points \c a, \c b
 * \c c.
 */
template <typename T>
Point<3,T> naive_closest_point_to_triangle(const Point<3,T>& p,
                                           const Point<3,T>& a,
                                           const Point<3,T>& b,
                                           const Point<3,T>& c) {
  
  typedef Point<3,T> point_type;

  // obtain plane of the triangle
  Plane pi(a,b,c);

  // get point in the plane closest to p
  point_type q = closest_ponint_to_plane(p,pi);
  
  // return if point is within the triangle
  if (is_point_in_triangle2(q, a, b, c))
    return q;
  
  // else get the closest point taking into account all edges
  
  // first edge
  q = closest_point_to_segment(p, a, b);
  T d = (q-p).sq_norm();
  
  // second edge
  point_type r = closest_point_to_segment(p, b, c);

  T d2 = (r-p).sq_norm();
  if (d2 < d) {
    q = r;
    d = d2;
  }
  
  // third edge
  r = closest_point_to_segment(p,c,a);
  d2 = (r-p).sq_norm();
  if (d2 < d)
    q = r;
  
  // return closest point
  return q;
}

__END_AKANTU__


#endif /* __AKANTU_GEOMETRY_HH__ */
