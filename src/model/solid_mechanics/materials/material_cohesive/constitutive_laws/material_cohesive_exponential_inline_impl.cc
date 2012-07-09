/**
 * @file   material_exponential_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:57:43 2010
 *
 * @brief  Implementation of the inline functions of the material exponential
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
// template<UInt dim>
// void  MaterialCohesiveExponential::computeTangentStiffness( Vector<Real> & normal, Real * tangent) {

// AKANTU_DEBUG_IN();

   // Vector<Real>::const_iterator<types::Matrix> normal_it =
   //    normal.begin(spatial_dimension,1);

   // for() {
   //   types::Matrix nn(spatial_dimension, spatial_dimension);
   //   nn.mul<false, true>(*normal_it, *normal_it);

   //   types::Matrix I(spatial_dimension, spatial_dimension);
   //   I.eye(beta*beta);
   //   nn *= (1-beta*beta);
   //   I += nn;


 
     

   // }


  // UInt n = (dim * (dim - 1) / 2 + dim);

  // Real Ep = E/((1+nu)*(1-2*nu));
  // Real Miiii = Ep * (1-nu);
  // Real Miijj = Ep * nu;
  // Real Mijij = Ep * (1-2*nu) * .5;

  // tangent[0 * n + 0] = Miiii;

  // // test of dimension should by optimized out by the compiler due to the template
  // if(dim >= 2) {
  //   tangent[1 * n + 1] = Miiii;
  //   tangent[0 * n + 1] = Miijj;
  //   tangent[1 * n + 0] = Miijj;

  //   tangent[(n - 1) * n + (n - 1)] = Mijij;
  // }

  // if(dim == 3) {
  //   tangent[2 * n + 2] = Miiii;
  //   tangent[0 * n + 2] = Miijj;
  //   tangent[1 * n + 2] = Miijj;
  //   tangent[2 * n + 0] = Miijj;
  //   tangent[2 * n + 1] = Miijj;

  //   tangent[3 * n + 3] = Mijij;
  //   tangent[4 * n + 4] = Mijij;
  // }

//  AKANTU_DEBUG_OUT();

// }

