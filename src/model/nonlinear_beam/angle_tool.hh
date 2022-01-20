#ifndef __ANGLE_TOOL_HPP__
#define __ANGLE_TOOL_HPP__

#include "aka_types.hh"
#include "aka_array.hh"
#include <cmath>
#include <iostream>

namespace akantu {


inline auto compute_rotation_angle(const Matrix<Real> &R) {
  //! From a rotation matrix, computes the rotation angle

  auto cos_val = (R.trace() - 1) / 2;
  if (cos_val < -1) {
    cos_val = -1;
  }
  if (cos_val > 1) {
    cos_val = 1;
  }

  auto theta = std::acos(cos_val);
  if (std::isnan(theta)) {
    std::cerr << R << std::endl;
    throw std::runtime_error(
        "in angle calculation theta is not a number: abort");
    // R, (R.trace() - 1) / 2);
  }
  return theta;
}







}

#endif //__ANGLE_TOOL_HPP__
