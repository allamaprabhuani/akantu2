#include "angle_tool.hh"

Vector3d set_axis(int i) {
  Vector3d res = Vector3d::Zero();
  res[i] = 1;
  return res;
}

Vector3d OO;
Vector3d e1 = set_axis(0);
Vector3d e2 = set_axis(1);
Vector3d e3 = set_axis(2);
