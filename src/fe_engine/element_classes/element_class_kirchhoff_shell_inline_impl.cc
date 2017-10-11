/**
 * @file   element_class_kirchhoff_shell_inline_impl.cc
 *
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation: Fri Jul 04 2014
 * @date last modification: Sun Oct 19 2014
 *
 * @brief  Element class Kirchhoff Shell
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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
#include "element_class_structural.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_CLASS_KIRCHHOFF_SHELL_INLINE_IMPL_CC__
#define __AKANTU_ELEMENT_CLASS_KIRCHHOFF_SHELL_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_STRUCTURAL_INTERPOLATION_TYPE_PROPERTY(_itp_kirchhoff_shell,
                                                     _itp_lagrange_triangle_3, 3, 6);
AKANTU_DEFINE_STRUCTURAL_ELEMENT_CLASS_PROPERTY(_kirchhoff_shell,
                                                _gt_triangle_3,
                                                _itp_kirchhoff_shell,
                                                _triangle_3, _ek_structural, 3,
                                                _git_triangle, 2);

/* -------------------------------------------------------------------------- */
// cf. element_class_bernoulli...
// element_class_triangle_3_inline... (pour compute Jacobian)
template <>
inline void InterpolationElement<_itp_kirchhoff_shell>::computeShapes(
    const Vector<Real> & /*natural_coords*/, Matrix<Real> & /*N*/,
    const Matrix<Real> & /*projected_coord*/) {
//   // projected_coord (x1 x2 x3) (y1 y2 y3)

//   // natural coordinate
//   Real xi = natural_coords(0);
//   Real eta = natural_coords(1);

//   Real x21 = projected_coord(0, 0) - projected_coord(0, 1); // x1-x2
//   Real x32 = projected_coord(0, 1) - projected_coord(0, 2);
//   Real x13 = projected_coord(0, 2) - projected_coord(0, 0);

//   Real y21 = projected_coord(1, 0) - projected_coord(1, 1); // y1-y2
//   Real y32 = projected_coord(1, 1) - projected_coord(1, 2);
//   Real y13 = projected_coord(1, 2) - projected_coord(1, 0);

//   /* Real x21=projected_coord(0,1)-projected_coord(0,0);
//   Real x32=projected_coord(0,2)-projected_coord(0,1);
//   Real x13=projected_coord(0,0)-projected_coord(0,2);

//   Real y21=projected_coord(1,1)-projected_coord(1,0);
//   Real y32=projected_coord(1,2)-projected_coord(1,1);
//   Real y13=projected_coord(1,0)-projected_coord(1,2);*/

//   // natural triangle side length
//   Real L4 = sqrt(x21 * x21 + y21 * y21);
//   Real L5 = sqrt(x32 * x32 + y32 * y32);
//   Real L6 = sqrt(x13 * x13 + y13 * y13);

//   // sinus and cosinus
//   Real C4 = x21 / L4; // 1
//   Real C5 = x32 / L5; //-1/sqrt(2);
//   Real C6 = x13 / L6; // 0
//   Real S4 = y21 / L4; // 0;
//   Real S5 = y32 / L5; // 1/sqrt(2);
//   Real S6 = y13 / L6; //-1;

//   Real N1 = 1 - xi - eta;
//   Real N2 = xi;
//   Real N3 = eta;

//   Real P4 = 4 * xi * (1 - xi - eta);
//   Real P5 = 4 * xi * eta;
//   Real P6 = 4 * eta * (1 - xi - eta);

//   // switch (id) {
//   // case 0: { // N
//   //   N(0) = N1;
//   //   N(1) = N2;
//   //   N(2) = N3;
//   //   break;
//   // }
//   // case 1: { // Nwi2
//   //   N(0) = -(1 / 8) * P4 * L4 * C4 + (1 / 8) * P6 * L6 * C6;
//   //   N(1) = -(1 / 8) * P5 * L5 * C5 + (1 / 8) * P4 * L4 * C4;
//   //   N(2) = -(1 / 8) * P6 * L6 * C6 + (1 / 8) * P5 * L5 * C5;
//   //   break;
//   // }
//   // case 2: { // Nwi3
//   //   N(0) = -(1 / 8) * P4 * L4 * S4 + (1 / 8) * P6 * L6 * S6;
//   //   N(1) = -(1 / 8) * P5 * L5 * S5 + (1 / 8) * P4 * L4 * S4;
//   //   N(2) = -(1 / 8) * P6 * L6 * S6 + (1 / 8) * P5 * L5 * S5;
//   //   break;
//   // }
//   // case 3: { // Nxi1
//   //   N(0) = 3 / (2 * L4) * P4 * C4 - 3 / (2 * L6) * P6 * C6;
//   //   N(1) = 3 / (2 * L5) * P5 * C5 - 3 / (2 * L4) * P4 * C4;
//   //   N(2) = 3 / (2 * L6) * P6 * C6 - 3 / (2 * L5) * P5 * C5;
//   //   break;
//   // }
//   // case 4: { // Nxi2
//   //   N(0) = N1 - (3. / 4.) * P4 * C4 * C4 - (3. / 4.) * P6 * C6 * C6;
//   //   N(1) = N2 - (3. / 4.) * P5 * C5 * C5 - (3. / 4.) * P4 * C4 * C4;
//   //   N(2) = N3 - (3. / 4.) * P6 * C6 * C6 - (3. / 4.) * P5 * C5 * C5;
//   //   break;
//   // }
//   // case 5: { // Nxi3
//   //   N(0) = -(3. / 4.) * P4 * C4 * S4 - (3. / 4.) * P6 * C6 * S6;
//   //   N(1) = -(3. / 4.) * P5 * C5 * S5 - (3. / 4.) * P4 * C4 * S4;
//   //   N(2) = -(3. / 4.) * P6 * C6 * S6 - (3. / 4.) * P5 * C5 * S5;
//   //   break;
//   // }
//   // case 6: { // Nyi1
//   //   N(0) = 3 / (2 * L4) * P4 * S4 - 3 / (2 * L6) * P6 * S6;
//   //   N(1) = 3 / (2 * L5) * P5 * S5 - 3 / (2 * L4) * P4 * S4;
//   //   N(2) = 3 / (2 * L6) * P6 * S6 - 3 / (2 * L5) * P5 * S5;
//   //   break;
//   // }
//   // case 7: { // Nyi2
//   //   N(0) = -(3. / 4.) * P4 * C4 * S4 - (3. / 4.) * P6 * C6 * S6;
//   //   N(1) = -(3. / 4.) * P5 * C5 * S5 - (3. / 4.) * P4 * C4 * S4;
//   //   N(2) = -(3. / 4.) * P6 * C6 * S6 - (3. / 4.) * P5 * C5 * S5;
//   //   break;
//   // }
//   // case 8: { // Nyi3
//   //   N(0) = N1 - (3. / 4.) * P4 * S4 * S4 - (3. / 4.) * P6 * S6 * S6;
//   //   N(1) = N2 - (3. / 4.) * P5 * S5 * S5 - (3. / 4.) * P4 * S4 * S4;
//   //   N(2) = N3 - (3. / 4.) * P6 * S6 * S6 - (3. / 4.) * P5 * S5 * S5;
//   //   break;
//   // }
//   // }
}

/* -------------------------------------------------------------------------- */
template <>
inline void
InterpolationElement<_itp_kirchhoff_shell, _itk_structural>::computeDNDS(
    const Vector<Real> & natural_coords, Matrix<Real> & B,
    const Matrix<Real> & projected_coord) {

  // natural coordinate
  Real xi = natural_coords(0);
  Real eta = natural_coords(1);

  // projected_coord (x1 x2 x3) (y1 y2 y3)

  // donne juste pour pour le patch test 4_5_5 mais donne quelque changement de
  // signe dans la matrice de rotation

  Real x21 = projected_coord(0, 0) - projected_coord(0, 1); // x1-x2
  Real x32 = projected_coord(0, 1) - projected_coord(0, 2);
  Real x13 = projected_coord(0, 2) - projected_coord(0, 0);

  Real y21 = projected_coord(1, 0) - projected_coord(1, 1); // y1-y2
  Real y32 = projected_coord(1, 1) - projected_coord(1, 2);
  Real y13 = projected_coord(1, 2) - projected_coord(1, 0);

  // donne juste pour la matrice de rigidité... mais pas pour le patch test
  // 4_5_5

  /* Real x21=projected_coord(0,1)-projected_coord(0,0);
    Real x32=projected_coord(0,2)-projected_coord(0,1);
    Real x13=projected_coord(0,0)-projected_coord(0,2);

    Real y21=projected_coord(1,1)-projected_coord(1,0);
    Real y32=projected_coord(1,2)-projected_coord(1,1);
    Real y13=projected_coord(1,0)-projected_coord(1,2);*/

  // natural triangle side length
  Real L4 = sqrt(x21 * x21 + y21 * y21);
  Real L5 = sqrt(x32 * x32 + y32 * y32);
  Real L6 = sqrt(x13 * x13 + y13 * y13);

  // sinus and cosinus
  Real C4 = x21 / L4;
  Real C5 = x32 / L5;
  Real C6 = x13 / L6;
  Real S4 = y21 / L4;
  Real S5 = y32 / L5;
  Real S6 = y13 / L6;

  Real dN1xi = -1;
  Real dN2xi = 1;
  Real dN3xi = 0;

  Real dN1eta = -1;
  Real dN2eta = 0;
  Real dN3eta = 1;

  Real dP4xi = 4 - 8 * xi - 4 * eta;
  Real dP5xi = 4 * eta;
  Real dP6xi = -4 * eta;

  Real dP4eta = -4 * xi;
  Real dP5eta = 4 * xi;
  Real dP6eta = 4 - 4 * xi - 8 * eta;

  // N'xi
  auto Np00 = dN1xi;
  auto Np01 = dN2xi;
  auto Np02 = dN3xi;
  //   N'eta
  auto Np10 = dN1eta;
  auto Np11 = dN2eta;
  auto Np12 = dN3eta;

  // Nxi1'xi
  auto Nx1p00 = 3. / (2 * L4) * dP4xi * C4 - 3. / (2. * L6) * dP6xi * C6;
  auto Nx1p01 = 3. / (2 * L5) * dP5xi * C5 - 3. / (2. * L4) * dP4xi * C4;
  auto Nx1p02 = 3. / (2 * L6) * dP6xi * C6 - 3. / (2. * L5) * dP5xi * C5;
  //    Nxi1'eta
  auto Nx1p10 = 3. / (2 * L4) * dP4eta * C4 - 3. / (2. * L6) * dP6eta * C6;
  auto Nx1p11 = 3. / (2 * L5) * dP5eta * C5 - 3. / (2. * L4) * dP4eta * C4;
  auto Nx1p12 = 3. / (2 * L6) * dP6eta * C6 - 3. / (2. * L5) * dP5eta * C5;

  // Nxi2'xi
  auto Nx2p00 = -1 - (3. / 4.) * dP4xi * C4 * C4 - (3. / 4.) * dP6xi * C6 * C6;
  auto Nx2p01 = 1 - (3. / 4.) * dP5xi * C5 * C5 - (3. / 4.) * dP4xi * C4 * C4;
  auto Nx2p02 = -(3. / 4.) * dP6xi * C6 * C6 - (3. / 4.) * dP5xi * C5 * C5;
  //    Nxi2'eta
  auto Nx2p10 =
      -1 - (3. / 4.) * dP4eta * C4 * C4 - (3. / 4.) * dP6eta * C6 * C6;
  auto Nx2p11 = -(3. / 4.) * dP5eta * C5 * C5 - (3. / 4.) * dP4eta * C4 * C4;
  auto Nx2p12 = 1 - (3. / 4.) * dP6eta * C6 * C6 - (3. / 4.) * dP5eta * C5 * C5;

  // Nxi3'xi
  auto Nx3p00 = -(3. / 4.) * dP4xi * C4 * S4 - (3. / 4.) * dP6xi * C6 * S6;
  auto Nx3p01 = -(3. / 4.) * dP5xi * C5 * S5 - (3. / 4.) * dP4xi * C4 * S4;
  auto Nx3p02 = -(3. / 4.) * dP6xi * C6 * S6 - (3. / 4.) * dP5xi * C5 * S5;
  //  Nxi3'eta
  auto Nx3p10 = -(3. / 4.) * dP4eta * C4 * S4 - (3. / 4.) * dP6eta * C6 * S6;
  auto Nx3p11 = -(3. / 4.) * dP5eta * C5 * S5 - (3. / 4.) * dP4eta * C4 * S4;
  auto Nx3p12 = -(3. / 4.) * dP6eta * C6 * S6 - (3. / 4.) * dP5eta * C5 * S5;

  // Nyi1'xi
  auto Ny1p00 = 3 / (2 * L4) * dP4xi * S4 - 3 / (2 * L6) * dP6xi * S6;
  auto Ny1p01 = 3 / (2 * L5) * dP5xi * S5 - 3 / (2 * L4) * dP4xi * S4;
  auto Ny1p02 = 3 / (2 * L6) * dP6xi * S6 - 3 / (2 * L5) * dP5xi * S5;
  //    Nyi1'eta
  auto Ny1p10 = 3 / (2 * L4) * dP4eta * S4 - 3 / (2 * L6) * dP6eta * S6;
  auto Ny1p11 = 3 / (2 * L5) * dP5eta * S5 - 3 / (2 * L4) * dP4eta * S4;
  auto Ny1p12 = 3 / (2 * L6) * dP6eta * S6 - 3 / (2 * L5) * dP5eta * S5;

  // Nyi2'xi
  auto Ny2p00 = -(3. / 4.) * dP4xi * C4 * S4 - (3. / 4.) * dP6xi * C6 * S6;
  auto Ny2p01 = -(3. / 4.) * dP5xi * C5 * S5 - (3. / 4.) * dP4xi * C4 * S4;
  auto Ny2p02 = -(3. / 4.) * dP6xi * C6 * S6 - (3. / 4.) * dP5xi * C5 * S5;
  //  Nyi2'eta
  auto Ny2p10 = -(3. / 4.) * dP4eta * C4 * S4 - (3. / 4.) * dP6eta * C6 * S6;
  auto Ny2p11 = -(3. / 4.) * dP5eta * C5 * S5 - (3. / 4.) * dP4eta * C4 * S4;
  auto Ny2p12 = -(3. / 4.) * dP6eta * C6 * S6 - (3. / 4.) * dP5eta * C5 * S5;

  // Nyi3'xi
  auto Ny3p00 =
      dN1xi - (3. / 4.) * dP4xi * S4 * S4 - (3. / 4.) * dP6xi * S6 * S6;
  auto Ny3p01 =
      dN2xi - (3. / 4.) * dP5xi * S5 * S5 - (3. / 4.) * dP4xi * S4 * S4;
  auto Ny3p02 =
      dN3xi - (3. / 4.) * dP6xi * S6 * S6 - (3. / 4.) * dP5xi * S5 * S5;
  //      Nyi3'eta
  auto Ny3p10 =
      dN1eta - (3. / 4.) * dP4eta * S4 * S4 - (3. / 4.) * dP6eta * S6 * S6;
  auto Ny3p11 =
      dN2eta - (3. / 4.) * dP5eta * S5 * S5 - (3. / 4.) * dP4eta * S4 * S4;
  auto Ny3p12 =
      dN3eta - (3. / 4.) * dP6eta * S6 * S6 - (3. / 4.) * dP5eta * S5 * S5;

  // clang-format off
  //       0    1                2                3                4       5     6     7                8                9               10   11    12    13               14               15               16
  B = {{Np00,   0.,              0.,              0.,              0.,     0., Np01,   0.,              0.,              0.,              0.,  0., Np02,   0.,              0.,              0.,              0., 0.},  // 0
       {  0., Np10,              0.,              0.,              0.,     0.,   0., Np11,              0.,              0.,              0.,  0.,   0., Np12,              0.,              0.,              0., 0.},  // 1
       {Np10, Np00,              0.,              0.,              0.,     0., Np11, Np01,              0.,              0.,              0.,  0., Np12, Np02,              0.,              0.,              0., 0.},  // 2
       {  0.,   0.,          Nx1p00,          Nx2p00,          Nx3p00,     0.,   0.,   0.,          Nx1p01,          Nx2p01,          Nx3p01,  0.,   0.,   0.,          Nx1p02,          Nx2p02,          Nx3p02, 0.},  // 3
       {  0.,   0.,          Ny1p10,          Ny2p10,          Ny3p10,     0.,   0.,   0.,          Ny1p11,          Ny2p11,          Ny3p11,  0.,   0.,   0.,          Ny1p12,          Ny2p12,          Ny3p12, 0.},  // 4
       {  0.,   0., Nx1p10 + Ny1p00, Nx2p10 + Ny2p00, Nx3p10 + Ny3p00,     0.,   0.,   0., Nx1p11 + Ny1p01, Nx2p11 + Ny2p01, Nx3p11 + Ny3p01,  0.,   0.,   0., Nx1p12 + Ny1p02, Nx2p12 + Ny2p02, Nx3p12 + Ny3p02, 0.}};  // 5
  // clang-format on
}


/* -------------------------------------------------------------------------- */
template <>
inline void ElementClass<_kirchhoff_shell,_ek_structural>::computeShapeDerivatives(
    const Matrix<Real> & /*natural_coord*/, Tensor3<Real> & /*shape_deriv*/,
    const Matrix<Real> & /*real_nodal_coord*/) {

  /// TO BE CONTINUED and moved in a _tmpl.hh
  //UInt spatial_dimension = real_nodal_coord.cols();
//   UInt nb_nodes = real_nodal_coord.rows();

//   const UInt projected_dim = natural_coord.rows();
//   Matrix<Real> rotation_matrix(real_nodal_coord);
//   Matrix<Real> rotated_nodal_coord(real_nodal_coord);
//   Matrix<Real> projected_nodal_coord(natural_coord);
//   /* ------------------------------------------------------------------------ */
//   Matrix<Real> Pe(real_nodal_coord);
//   Matrix<Real> Pg(real_nodal_coord);
//   Matrix<Real> inv_Pg(real_nodal_coord);

//   /// compute matrix Pe
//   Pe.eye();

//   // /// compute matrix Pg
//   // Pg(0) = real_nodal_coord(1) - real_nodal_coord(0);
//   // Pg(1) = real_nodal_coord(2) - real_nodal_coord(0);
//   // Pg(2).crossProduct(Pg(0), Pg(1));

//   // /// compute inverse of Pg
//   // inv_Pg.inverse(Pg);

//   /// compute rotation matrix
//   // rotation_matrix=Pe*inv_Pg;
//   rotation_matrix.eye();

//   /* ------------------------------------------------------------------------ */
//   rotated_nodal_coord.mul<false, false>(rotation_matrix, real_nodal_coord);

//   UInt nb_points = shape_deriv.size(2);

//   for (UInt i = 0; i < projected_dim; ++i) {
//     for (UInt j = 0; j < nb_points; ++j) {
//       projected_nodal_coord(i, j) = rotated_nodal_coord(i, j);
//     }
//   }

//   Tensor3<Real> dnds(projected_dim, nb_nodes, natural_coord.cols());
//   Tensor3<Real> J(projected_dim, projected_dim, natural_coord.cols());

//   parent_element::computeDNDS(natural_coord, dnds);
//   parent_element::computeJMat(dnds, projected_nodal_coord, J);

//   for (UInt p = 0; p < nb_points; ++p) {
//     Matrix<Real> shape_deriv_p = shape_deriv(p);

//     interpolation_element::computeDNDS(natural_coord(p), shape_deriv_p,
//                                        projected_nodal_coord);

//     Matrix<Real> dNdS = shape_deriv_p;
//     Matrix<Real> inv_J(projected_dim, projected_dim);
//     inv_J.inverse(J(p));
//     shape_deriv_p.mul<false, false>(inv_J, dNdS);
//   }
}

}  // akantu

#endif /* __AKANTU_ELEMENT_CLASS_KIRCHHOFF_SHELL_INLINE_IMPL_CC__ */
