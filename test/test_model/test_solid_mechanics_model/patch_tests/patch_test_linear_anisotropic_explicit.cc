/**
 * @file   patch_test_explicit_anisotropic.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Sat Apr 16 2011
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  patch test for elastic material in solid mechanics model
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "patch_test_linear_fixture.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

// Stiffness tensor, rotated by hand
Real C[3][3][3]
      [3] = {{{{112.93753505, 1.85842452538e-10, -4.47654358027e-10},
               {1.85847317471e-10, 54.2334345331, -3.69840984824},
               {-4.4764768395e-10, -3.69840984824, 56.848605217}},
              {{1.85847781609e-10, 25.429294233, -3.69840984816},
               {25.429294233, 3.31613847493e-10, -8.38797920011e-11},
               {-3.69840984816, -8.38804581349e-11, -1.97875715813e-10}},
              {{-4.47654358027e-10, -3.69840984816, 28.044464917},
               {-3.69840984816, 2.09374961813e-10, 9.4857455224e-12},
               {28.044464917, 9.48308098714e-12, -2.1367885239e-10}}},
             {{{1.85847781609e-10, 25.429294233, -3.69840984816},
               {25.429294233, 3.31613847493e-10, -8.38793479119e-11},
               {-3.69840984816, -8.38795699565e-11, -1.97876381947e-10}},
              {{54.2334345331, 3.31617400207e-10, 2.09372075233e-10},
               {3.3161562385e-10, 115.552705733, -3.15093728886e-10},
               {2.09372075233e-10, -3.15090176173e-10, 54.2334345333}},
              {{-3.69840984824, -8.38795699565e-11, 9.48219280872e-12},
               {-8.38795699565e-11, -3.1509195253e-10, 25.4292942335},
               {9.48441325477e-12, 25.4292942335, 3.69840984851}}},
             {{{-4.47653469848e-10, -3.69840984816, 28.044464917},
               {-3.69840984816, 2.09374073634e-10, 9.48752187924e-12},
               {28.044464917, 9.48552347779e-12, -2.1367885239e-10}},
              {{-3.69840984824, -8.3884899027e-11, 9.48219280872e-12},
               {-8.3884899027e-11, -3.150972816e-10, 25.4292942335},
               {9.48041645188e-12, 25.4292942335, 3.69840984851}},
              {{56.848605217, -1.97875493768e-10, -2.13681516925e-10},
               {-1.97877270125e-10, 54.2334345333, 3.69840984851},
               {-2.13683293282e-10, 3.69840984851, 112.93753505}}}};

/* -------------------------------------------------------------------------- */
TYPED_TEST(TestPatchTestLinear, AnisotropicExplicit) {
  this->initModel(_explicit_lumped_mass, "material_anisotropic.dat");

  const auto & coordinates = this->mesh->getNodes();
  auto & displacement = this->model->getDisplacement();

  // set the position of all nodes to the static solution
  for (auto && tuple : zip(make_view(coordinates, this->dim),
                           make_view(displacement, this->dim))) {
    this->setDisplacement(std::get<1>(tuple), std::get<0>(tuple));
  }

  for (UInt s = 0; s < 100; ++s) {
    this->model->solveStep();
  }

  auto ekin = this->model->getEnergy("kinetic");
  EXPECT_NEAR(0, ekin, 1e-16);

  this->checkDisplacements();
  this->checkStrains();
  this->checkStresses([&](const Matrix<Real> & pstrain) {
    auto strain = (pstrain + pstrain.transpose()) / 2.;
    decltype(strain) stress(this->dim, this->dim);

    for (UInt i = 0; i < this->dim; ++i) {
      for (UInt j = 0; j < this->dim; ++j) {
        stress(i, j) = 0;
        for (UInt k = 0; k < this->dim; ++k) {
          for (UInt l = 0; l < this->dim; ++l) {
            stress(i, j) += C[i][j][k][l] * strain(k, l);
          }
        }
      }
    }
    return stress;
  });
}
