/**
 * @file   test_solid_mechanics_model_implicit_2d.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Mar 03 16:09:49 2011
 *
 * @brief  test of traction in implicit
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

/* -------------------------------------------------------------------------- */
#include <limits>
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "static_communicator.hh"
#include "mesh_partition_scotch.hh"

/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_SCOTCH
#include "mesh_partition_scotch.hh"
#endif

#define bar_length 1
#define bar_height 1


/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  akantu::debug::setDebugLevel(akantu::dblWarning);
  akantu::initialize(argc, argv);

  akantu::UInt spatial_dimension = 2;

  akantu::Mesh mesh(spatial_dimension);

  akantu::StaticCommunicator & comm = akantu::StaticCommunicator::getStaticCommunicator();
  akantu::Int psize = comm.getNbProc();
  akantu::Int prank = comm.whoAmI();

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    akantu::MeshIOMSH mesh_io;
    mesh_io.read("square_implicit2.msh", mesh);

    partition = new akantu::MeshPartitionScotch(mesh, spatial_dimension);
    //   partition->reorder();
    partition->partitionate(psize);
  }
  akantu::SolidMechanicsModel model(mesh);

  /// model initialization
  model.initParallel(partition);
  delete partition;

  model.initFull("material.dat", akantu::_static);

  if (prank == 0)
    std::cout << model.getMaterial(0) << std::endl;

  /// boundary conditions
  const  akantu::Array<akantu::Real> & position = mesh.getNodes();
  akantu::Array<bool> & boundary = model.getBoundary();
  akantu::Array<akantu::Real> & displacment = model.getDisplacement();

  akantu::UInt nb_nodes = model.getFEM().getMesh().getNbNodes();
  for (akantu::UInt n = 0; n < nb_nodes; ++n) {
    if(position(n,0) < akantu::Math::getTolerance()) boundary(n,0) = true;
    if(position(n,1) < akantu::Math::getTolerance()) boundary(n,1) = true;

    if(std::abs(position(n,0) - bar_length) < akantu::Math::getTolerance()) {
      boundary(n,0) = true;
      displacment(n,0) = 0.1;
    }
  }

  akantu::UInt count = 0;
  model.updateResidual();


  model.setBaseName("implicit_2d");
  model.addDumpField("displacement");
  model.addDumpField("mass"        );
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("force"       );
  model.addDumpField("residual"    );
  model.addDumpField("stress"      );
  model.addDumpField("strain"      );

  model.computeStresses();
  model.dump();

  akantu::Real norm;
  model.assembleStiffnessMatrix();
  while(!model.testConvergenceResidual(1e-3, norm) && (count < 100)) {
    if (prank == 0)
      std::cout << "Iter : " << ++count << " - residual norm : " << norm << std::endl;

    model.solveStatic();
    model.updateResidual();
  }

  model.computeStresses();
  model.dump();

  akantu::finalize();

  return EXIT_SUCCESS;
}
