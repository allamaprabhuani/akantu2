/**
 * @file   test_gradient_bernoulli_beam_2.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Tue Apr  5 17:19:48 2011
 *
 * @brief Test of the gradient on the type _bernoulli_beam_2
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
#include <cstdlib>
#include <fstream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "fem.hh"
#include "integrator_gauss.hh"
#include "shape_linked.hh"

/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]){

  Vector<Real> displ_on_nodes;
  Vector<Real> grad_on_quad;
  UInt k=0;

  
  displ_on_nodes=(2,0,0,4,0,0,6,0,0,8,0,0); // The displacement fields to interpolate
  UInt size =displ_on_nodes.getSize();

 // The gradient is only defined for du/dx and -d(dv)/dx^2 used in the definition of N and M stresses.
   

  while (k<size) {
    
    gradientOnControlPoints(displ_on_nodes,grad_on_quad,0,k,k,
			       size,false,_bernoulli_beam_2);
    ++k;

    gradientOnControlPoints(displ_on_nodes, grad_on_quad,1,(k-1),k,
				size,false,_bernoulli_beam_2);

    gradientOnControlPoints(displ_on_nodes, grad_on_quad,2,k,k,
				size,true,_bernoulli_beam_2);
 
  }

 std::ofstream my_file("out.txt");
  my_file << displ_on_nodes << std::endl;
  my_file << grad_on_quad << std::endl;

  return EXIT_SUCCESS;

} 
