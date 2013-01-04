/**
 * @file   aka_sphere.hh
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @date   Mon Jun 18 10:20:00 2012
 *
 * @brief  bounding sphere classes
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

#ifndef __AKANTU_SPHERE_HH__
#define __AKANTU_SPHERE_HH__

#include <iostream>

#include "aka_common.hh"
#include "aka_point.hh"
#include "aka_bounding_box.hh"

__BEGIN_AKANTU__

using std::cout;
using std::endl;

// forward declarations
class Sphere;
//void update_sphere(Sphere&, const typename Sphere::point_type&);


//! Sphere class
class Sphere : public Bounding_volume<3> {
  
public:
  
  static constexpr int dim_ = 3;
  
  typedef Bounding_volume<dim_> base_type;
  typedef typename base_type::point_type point_type;

  typedef typename point_type::value_type value_type;
  typedef BoundingBox<dim_> aabb_type;


  static int dim()
  { return dim_; }
  
  Sphere(const point_type& c = point_type(), value_type r = value_type()) : base_type(), c_(c), r_(r) {}
  
  virtual base_type* combine(const base_type& b) const {
    
    const Sphere* sp = dynamic_cast<const Sphere*>(&b);
    assert(sp != nullptr);
    const Sphere& s0 = *sp;
    Sphere r(s0);
    r += *this;
    return new Sphere(r);
  }
  
  virtual std::ostream& print(std::ostream& os) const {
    os<<"Sphere["<<c_<<", "<<r_<<"]";
    return os;
  }
  
  aabb_type aabb() const {

    point_type o = r_*point_type(1.);
    return aabb_type(c_ - o, c_ + o);
  }
  
  
  point_type const& center() const
  { return c_; }
  
  value_type const& radius() const
  { return r_; }
  
  //! Use in generic code as comparative measure of how big the sphere is
  value_type measure() const
  { return pow(r_,3); }
  
  //! Sphere intersection test
  bool operator&(const Sphere& s) {
    point_type d = c_ - s.c_;
    value_type sq_dist = d.sq_norm();
    return sq_dist <= pow(r_ + s.r_,2);
  }

  //! Grow sphere if point lies outside of it
  Sphere& operator+=(const point_type& p) {
    
    point_type diff = p - c_;
    value_type sq_norm = diff.sq_norm();
    
    if (sq_norm > r_*r_) {
      
      value_type norm = sqrt(sq_norm);
      value_type new_r = 0.5*(r_ + norm);
      value_type scalar = (new_r - r_) / norm;
      r_ = new_r;
      c_ += scalar * diff;
    }
    return *this;
  }
  
  //! Determine the sphere that encloses both spheres
  Sphere& operator+=(const Sphere s) {
    
    point_type diff = s.c_ - c_;
    value_type sq_norm = diff.sq_norm();
    
    // one sphere is contained within the other
    if (pow(s.r_ - r_, 2) >= sq_norm) {
      if(s.r_ >= r_)
        this->operator=(s);
      // else do nothing, as the current sphere is bigger
      // and no further changes are required
    }
    // else spheres partially overlapping or disjoint
    else {
      
      // compute new radius
      value_type norm = sqrt(sq_norm);
      value_type tmp = r_;
      r_ = 0.5 * (norm + r_ + s.r_);
      if (norm > 2*std::numeric_limits<Real>::epsilon())
        c_ += ((r_ - tmp) / norm) * diff;
    }
    return *this;
  }
  
  
  bool operator&(const point_type& p) const
  { return (p - c_).sq_norm() < r_*r_; }
  
  bool operator&(const Sphere& s) const
  { return (c_ - s.c_).sq_norm() <= pow(r_ + s.r_,2.); }  
  
private:
  
  point_type c_;  //!< Sphere center
  Real r_;        //!< Sphere radius
};


Sphere operator+(const Sphere& s1, const Sphere& s2) {
  Sphere r(s1);
  return r += s2;
}


std::ostream& operator<<(std::ostream& os, const Sphere& s) {
  os<<"Sphere["<<s.center()<<", "<<s.radius()<<"]";
  return os;
}


// Ritter algorithm
// J. Ritter, Graphics gems, Academic Press Professional, Inc., San Diego, CA, USA, 1990, Ch. An efficient bounding sphere, pp. 301–303.
// URL http://dl.acm.org/citation.cfm?id=90767.90836

template <class point_container>
std::pair<size_t, size_t> extreme_points(const point_container& pts) {
  
  typedef typename point_container::value_type point_type;
  typedef typename point_type::value_type value_type;
  typedef typename point_container::const_iterator const_iterator;
  
  size_t min[] = { 0, 0, 0 };
  size_t max[] = { 0, 0, 0 };
  
  // loop over container points to find extremal points
  for (size_t i=1; i<pts.size(); ++i) {
    
    const point_type& p = pts[i];
    
    // loop over coordinates
    for (int j=0; j<point_type::dim(); ++j) {
      
      // check if new point is minimum
      if (p[j] < pts[min[j]][j])
        min[j] = i;
      // check if new point is maximum
      else if (p[j] > pts[max[j]][j])
        max[j] = i;
    }
  }
  
  // pick the pair of the longest distance
  size_t m=0, M=0;
  value_type sq_norm = value_type();
  
  for (int i=0; i<point_type::dim(); ++i) {
    
    point_type diff = pts[max[i]] - pts[min[i]];
    value_type new_sq_norm = diff.sq_norm();
    
    if (new_sq_norm > sq_norm) {
      m = min[i];
      M = max[i];
      sq_norm = new_sq_norm;
    }
  }
  
  return std::make_pair(m,M);
}


template <class point_container>
Sphere bounding_sphere(const point_container& pts) {
  
  assert(!pts.empty());
  
  typedef typename point_container::value_type point_type;
  typedef typename point_type::value_type value_type;
  
  // find extreme points on axis-aligned bounding box to construct
  // first approximation of the sphere
  std::pair<size_t, size_t> mM = extreme_points(pts);
  
  // compute center and radius of the sphere
  const point_type &m = pts[mM.first];
  const point_type &M = pts[mM.second];
  point_type c = 0.5*(m+M);
  value_type r = sqrt((M-c).sq_norm());
  
  // construct sphere
  Sphere s(c,r);
  
  // second pass: update the sphere so that all points lie inside
  for (size_t i=0; i<pts.size(); ++i)
    s += pts[i];
  
  return s;
}




__END_AKANTU__

#endif /* __AKANTU_SPHERE_HH__ */
