#include "angle_tool.hh"

Array<Real> set_axis(UInt i) {
  Array<Real> res(1,3)->zero();
  res[i] = 1;
  return res;
}

Array<Real> OO(1,3)->zero();
Array<Real> e1 = set_axis(0);
Array<Real> e2 = set_axis(1);
Array<Real> e3 = set_axis(2);
