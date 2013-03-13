/**
 * @file   structural_mechanics_model_boundary.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed May 25 15:21:50 2011
 *
 * @brief
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
#include "model.hh"
#include "structural_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::transferNMatrixToSymVoigtNMatrix<_bernoulli_beam_2>(Array<Real> & N_matrix) {
  AKANTU_DEBUG_IN();
  MyFEMType & fem = getFEMClass<MyFEMType>();
  UInt nb_nodes_per_element = getFEM().getMesh().getNbNodesPerElement(_bernoulli_beam_2);

  Array<Real>::const_iterator< Vector<Real> > shape_N0 = fem.getShapeFunctions().getShapes(_bernoulli_beam_2, _not_ghost, 0).begin(nb_nodes_per_element);
  Array<Real>::const_iterator< Vector<Real> > shape_M0 = fem.getShapeFunctions().getShapes(_bernoulli_beam_2, _not_ghost, 1).begin(nb_nodes_per_element);
  Array<Real>::const_iterator< Vector<Real> > shape_L0 = fem.getShapeFunctions().getShapes(_bernoulli_beam_2, _not_ghost, 2).begin(nb_nodes_per_element);
  Array<Real>::const_iterator< Vector<Real> > shape_Mp = fem.getShapeFunctions().getShapes(_bernoulli_beam_2, _not_ghost, 3).begin(nb_nodes_per_element);
  Array<Real>::const_iterator< Vector<Real> > shape_Lp = fem.getShapeFunctions().getShapes(_bernoulli_beam_2, _not_ghost, 4).begin(nb_nodes_per_element);

  N_matrix.clear();
  Array<Real>::iterator< Matrix<Real> > N_it = N_matrix.begin(nb_degree_of_freedom, nb_degree_of_freedom * nb_nodes_per_element);
  Array<Real>::iterator< Matrix<Real> > N_end = N_matrix.end(nb_degree_of_freedom, nb_degree_of_freedom * nb_nodes_per_element);

  for (;N_it != N_end; ++N_it, ++shape_N0, ++shape_M0, ++shape_L0, ++shape_Mp, ++shape_Lp) {
    Matrix<Real> & N = *N_it;
    const Vector<Real> & N0 = *shape_N0;
    const Vector<Real> & M0 = *shape_M0;
    const Vector<Real> & L0 = *shape_L0;
    const Vector<Real> & Mp = *shape_Mp;
    const Vector<Real> & Lp = *shape_Lp;

    N(0,0) = N0(0);
    N(0,3) = N0(1);

    N(1,1) = M0(0);
    N(1,2) = L0(0);
    N(1,4) = M0(1);
    N(1,5) = L0(1);

    N(2,1) = Mp(0);
    N(2,2) = Lp(0);
    N(2,4) = Mp(1);
    N(2,5) = Lp(1);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::transferNMatrixToSymVoigtNMatrix<_bernoulli_beam_3>(Array<Real> & N_matrix) {
  AKANTU_DEBUG_IN();

  ElementType type = _bernoulli_beam_3;

  MyFEMType & fem = getFEMClass<MyFEMType>();
  UInt nb_nodes_per_element = getFEM().getMesh().getNbNodesPerElement(type);

  Array<Real>::const_iterator< Vector<Real> > shape_N0 = fem.getShapeFunctions().getShapes(type, _not_ghost, 0).begin(nb_nodes_per_element);
  Array<Real>::const_iterator< Vector<Real> > shape_M0 = fem.getShapeFunctions().getShapes(type, _not_ghost, 1).begin(nb_nodes_per_element);
  Array<Real>::const_iterator< Vector<Real> > shape_L0 = fem.getShapeFunctions().getShapes(type, _not_ghost, 2).begin(nb_nodes_per_element);
  Array<Real>::const_iterator< Vector<Real> > shape_Mp = fem.getShapeFunctions().getShapes(type, _not_ghost, 3).begin(nb_nodes_per_element);
  Array<Real>::const_iterator< Vector<Real> > shape_Lp = fem.getShapeFunctions().getShapes(type, _not_ghost, 4).begin(nb_nodes_per_element);

  N_matrix.clear();
  Array<Real>::iterator< Matrix<Real> > N_it = N_matrix.begin(nb_degree_of_freedom, nb_degree_of_freedom * nb_nodes_per_element);
  Array<Real>::iterator< Matrix<Real> > N_end = N_matrix.end(nb_degree_of_freedom, nb_degree_of_freedom * nb_nodes_per_element);

  for (; N_it != N_end; ++N_it, ++shape_N0, ++shape_M0, ++shape_L0, ++shape_Mp, ++shape_Lp) {
    Matrix<Real> & N = *N_it;
    const Vector<Real> & N0 = *shape_N0;
    const Vector<Real> & M0 = *shape_M0;
    const Vector<Real> & L0 = *shape_L0;
    const Vector<Real> & Mp = *shape_Mp;
    const Vector<Real> & Lp = *shape_Lp;

    N(0,0)  =  N0(0);
    N(0,6)  =  N0(1);

    N(1,1)  =  M0(0);
    N(1,5)  =  L0(0);
    N(1,7)  =  M0(1);
    N(1,11) =  L0(1);

    N(2,2)  =  M0(0);
    N(2,4)  = -L0(0);
    N(2,8)  =  M0(1);
    N(2,10) = -L0(1);

    N(3,3)  =  N0(0);
    N(3,9)  =  N0(1);

    N(4,2)  =  Mp(0);
    N(4,4)  = -Lp(0);
    N(4,8)  =  Mp(1);
    N(4,10) = -Lp(1);

    N(5,1)  =  Mp(0);
    N(5,5)  =  Lp(0);
    N(5,7)  =  Mp(1);
    N(5,11) =  Lp(1);
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */

__END_AKANTU__
