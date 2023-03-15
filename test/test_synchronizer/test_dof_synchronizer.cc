/**
 * Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 * 
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 * 
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "communicator.hh"
#include "dof_synchronizer.hh"
#include "element_synchronizer.hh"
#include "mesh_io.hh"
#include "mesh_partition_scotch.hh"

/* -------------------------------------------------------------------------- */
#include "io_helper.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char * argv[]) {
  const Int spatial_dimension = 2;

  initialize(argc, argv);

  const auto & comm = akantu::Communicator::getStaticCommunicator();
  Int prank = comm.whoAmI();

  Mesh mesh(spatial_dimension);

  if (prank == 0)
    mesh.read("bar.msh");
  mesh.distribute();

  DOFManagerDefault dof_manager(mesh, "test_dof_manager");

  UInt nb_nodes = mesh.getNbNodes();

  /* ------------------------------------------------------------------------ */
  /* test the synchronization                                                 */
  /* ------------------------------------------------------------------------ */
  Array<Real> test_synchronize(nb_nodes, spatial_dimension, "Test vector");
  dof_manager.registerDOFs("test_synchronize", test_synchronize, _dst_nodal);

  const auto & equation_number =
      dof_manager.getLocalEquationsNumbers("test_synchronize");

  auto & dof_synchronizer = dof_manager.getSynchronizer();

  std::cout << "Synchronizing a dof vector" << std::endl;

  Array<Int> local_data_array(dof_manager.getLocalSystemSize(), 2);
  auto it_data = local_data_array.begin(2);

  for (Int local_dof = 0; local_dof < dof_manager.getLocalSystemSize();
       ++local_dof) {
    auto equ_number = equation_number(local_dof);
    Vector<Int> val;
    if (dof_manager.isLocalOrMasterDOF(equ_number)) {
      Int global_dof = dof_manager.localToGlobalEquationNumber(local_dof);
      val = {0, 1};
      val.array() += global_dof * 2;
    } else {
      val = {-1, -1};
    }

    it_data[local_dof] = val;
  }

  dof_synchronizer.synchronizeArray(local_data_array);

  auto test_data = [&]() -> void {
    auto it_data = local_data_array.begin(2);
    for (Int local_dof = 0; local_dof < dof_manager.getLocalSystemSize();
         ++local_dof) {
      Int equ_number = equation_number(local_dof);

      Vector<Int> exp_val;

      auto global_dof = dof_manager.localToGlobalEquationNumber(local_dof);

      if (dof_manager.isLocalOrMasterDOF(equ_number) ||
          dof_manager.isSlaveDOF(equ_number)) {
        exp_val = {0, 1};
        exp_val.array() += global_dof * 2;
      } else {
        exp_val = {-1, -1};
      }

      auto && val = it_data[local_dof];
      if (exp_val != val) {
        std::cerr << "Failed !" << prank << " DOF: " << global_dof << " - l"
                  << local_dof << " value:" << val << " expected: " << exp_val
                  << std::endl;
        exit(1);
      }
    }
  };

  test_data();

  if (prank == 0) {
    Array<Int> test_gathered(dof_manager.getSystemSize(), 2);
    dof_synchronizer.gather(local_data_array, test_gathered);

    local_data_array.set(-1);
    dof_synchronizer.scatter(local_data_array, test_gathered);
  } else {
    dof_synchronizer.gather(local_data_array);

    local_data_array.set(-1);
    dof_synchronizer.scatter(local_data_array);
  }

  test_data();

  finalize();

  return 0;
}
