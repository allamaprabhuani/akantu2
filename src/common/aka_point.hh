/**
 * @file   aka_point.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @author Alejandro Aragon <alejandro.aragon@epfl.ch>
 * @date   Thu Jan 12 15:43:51 2012
 *
 * @brief  class of a geometrical point
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

#ifndef __AKANTU_AKA_POINT_HH__
#define __AKANTU_AKA_POINT_HH__

#include "aka_common.hh"
#include <cmath>

#include <cassert>

__BEGIN_AKANTU__


const Real inf = std::numeric_limits<Real>::infinity();


template <int d, typename T = Real>
class Point {
  
public:
  
  typedef T value_type;
  
  static int dim()
  { return d; }
  
  Point() {
    for (UInt i=0; i<d; ++i)
      coord_[i] = value_type();
  }
  
  
  explicit Point(value_type const* coordinates) {
    for (UInt i=0; i<d; ++i)
      coord_[i] = coordinates[i];
  }
  
  // parameter constructor
  template <typename... Args>
  explicit Point(const Args&... args) {
    
    static_assert(sizeof...(Args) <= d , "*** ERROR *** Number of arguments exceeded point dension");
    
    std::fill_n(coord_, d, value_type());
    
    value_type coord[] = { args... };
    
    if (sizeof...(Args) != 0)
      for (size_t i=0; i<d; ++i)
        coord_[i] = i < sizeof...(Args) ? coord[i] : coord[sizeof...(Args) - 1];
  }
  
  /// less operator
  bool operator<(const Point & p) const {
    for (int i=0; i<d; ++i)
      if (coord_[i] < p[i])
        return true;
    return false;
  }
  
  /// bool equal operator
  bool operator==(const Point & p) const {
    for (int i=0; i<d; ++i)
      if (coord_[i] != p[i])
        return false;
    return true;
  }
  
  /// standard output stream operator
  friend std::ostream& operator <<(std::ostream &os, const Point &p) {
    os<<"{"<<p.coord_[0];
    for (int i=1; i<d; ++i)
      os<<","<<p.coord_[i];
    os<<"}";
    return os;
  }
  
public:
  
  //! Get copy of coordinate
  value_type operator[] (UInt index) const {
    assert(index < d);
    return coord_[index];
  }
  
  //! Get write access to coordinate
  value_type& operator[] (UInt index) {
    assert(index < d);
    return coord_[index];
  }
  
  Point& operator+=(const Point& p) {
    for (int i=0; i<d; ++i)
      coord_[i] += p.coord_[i];
    return *this;
  }
  
  Point& operator-=(const Point& p) {
    for (int i=0; i<d; ++i)
      coord_[i] -= p.coord_[i];
    return *this;
  }
  
  template <typename S>
  Point& operator*=(S s) {
    for (int i=0; i<d; ++i)
      coord_[i] *= s;
    return *this;
  }
  
  Point& normalize()
  { return (*this)*=(1/std::sqrt(sq_norm())); }
  
  value_type sq_norm() {
    value_type r = value_type();
    for (int i=0; i<d; ++i)
      r += std::pow(coord_[i],2);
    return r;
  }
  
private:
  
  Real coord_[d];
};


template <int d, typename T>
Point<d,T> operator+(const Point<d,T>& p, const Point<d,T>& q) {
  Point<d,T> r(p);
  return r += q;
}


template <int d, typename T>
Point<d,T> operator-(const Point<d,T>& p, const Point<d,T>& q) {
  Point<d,T> r(p);
  return r -= q;
}

// used operator* for dot product
template <int d, typename T>
typename Point<d,T>::value_type operator*(const Point<d,T>& p, const Point<d,T>& q) {
  
  typename Point<d,T>::value_type r = 0;
  for (int i=0; i<d; ++i)
    r += p[i] * q[i];
  return r;
}



template <int d, typename T>
Point<d,T> operator*(const Point<d,T>& p, typename Point<d,T>::value_type s) {
  Point<d,T> r(p);
  return r *= s;
}

template <int d, typename T>
Point<d,T> operator*(typename Point<d,T>::value_type s, const Point<d,T>& p) {
  Point<d,T> r(p);
  return r *= s;
}


// cross product
template <typename T>
Point<3,T> cross(const Point<3,T>& o, const Point<3,T>& p) {
  
  Point<3,T> r;
  for (int i=0; i<3; ++i)
    r[i] = o[(i+1)%3]*p[(i+2)%3] - o[(i+2)%3]*p[(i+1)%3];
  return r;
}


template <int d>
struct Bounding_volume {

  typedef Point<d> point_type;
  typedef typename point_type::value_type value_type;

  virtual value_type measure() const = 0;
  virtual std::ostream& print(std::ostream& os) const = 0;
  
  Real last_time_;
  point_type velocity_;
  point_type acceleration_;
  
  friend std::ostream& operator<<(std::ostream& os, const Bounding_volume& gv)
  { return gv.print(os); }
};




__END_AKANTU__

#endif /* __AKANTU_AKA_POINT_HH__ */
