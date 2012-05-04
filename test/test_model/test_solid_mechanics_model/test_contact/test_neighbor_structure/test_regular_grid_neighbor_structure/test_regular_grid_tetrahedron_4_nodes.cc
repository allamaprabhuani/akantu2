/**
 * @file   test_regular_grid_tetrahedron_4_nodes.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Oct 26 16:58:42 2010
 *
 * @brief  test regular grid neighbor structure for 3d case
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

#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "contact.hh"
#include "contact_rigid.hh"
#include "contact_neighbor_structure.hh"
#include "regular_grid_neighbor_structure.hh"

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"

#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);

  int dim = 3;

  /// load mesh
  Mesh my_mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("cubes.msh", my_mesh);

  /// build facet connectivity and surface id
  MeshUtils::buildFacets(my_mesh);
  MeshUtils::buildSurfaceID(my_mesh);
 
  unsigned int nb_nodes = my_mesh.getNbNodes();

  /// dump facet and surface information to paraview
#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperParaview dumper;
  dumper.SetMode(iohelper::TEXT);
  
  dumper.SetPoints(my_mesh.getNodes().values, dim, nb_nodes, "tetrahedron_4_nodes_test-surface-extraction");
  dumper.SetConnectivity((int*)my_mesh.getConnectivity(_tetrahedron_4).values,
   			 iohelper::TETRA1, my_mesh.getNbElement(_tetrahedron_4), iohelper::C_MODE);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  /// declaration of model
  SolidMechanicsModel  my_model(my_mesh);
  /// model initialization
  my_model.initVectors();
  // initialize the vectors
  memset(my_model.getForce().values,        0, 3*nb_nodes*sizeof(Real));
  memset(my_model.getVelocity().values,     0, 3*nb_nodes*sizeof(Real));
  memset(my_model.getAcceleration().values, 0, 3*nb_nodes*sizeof(Real));
  memset(my_model.getDisplacement().values, 0, 3*nb_nodes*sizeof(Real));

  my_model.initModel();  
  my_model.readMaterials("material.dat");
  my_model.initMaterials();

  Real time_step = my_model.getStableTimeStep();
  my_model.setTimeStep(time_step/10.);

  my_model.assembleMassLumped();

   /// contact declaration
  Contact * contact = Contact::newContact(my_model, 
					  _ct_rigid, 
					  _cst_expli, 
					  _cnst_regular_grid);

  ContactRigid * my_contact = dynamic_cast<ContactRigid *>(contact);

  my_contact->initContact(false);

  Surface master = 0;
  Surface impactor = 1;

  my_contact->addMasterSurface(master);
  my_contact->addImpactorSurfaceToMasterSurface(impactor, master);

  my_model.updateCurrentPosition(); // neighbor structure uses current position for init
  my_contact->initNeighborStructure(master);
  
  const NodesNeighborList & my_neighbor_list = dynamic_cast<const NodesNeighborList &>(my_contact->getContactSearch().getContactNeighborStructure(master).getNeighborList());

  UInt nb_nodes_neigh = my_neighbor_list.impactor_nodes.getSize();
  Vector<UInt> impact_nodes = my_neighbor_list.impactor_nodes;
  UInt * impact_nodes_val = impact_nodes.values;

  /// define output file for testing
  std::ofstream test_output;
  test_output.open("test_regular_grid_tetrahedron_4_nodes.out");

  /// print impactor nodes
  test_output << "we have " << nb_nodes_neigh << " impactor nodes:";
  for (UInt i = 0; i < nb_nodes_neigh; ++i) {
    test_output << " " << impact_nodes_val[i];
  }
  test_output << std::endl;

  UInt * master_nodes_offset_val = my_neighbor_list.master_nodes_offset.values;
  UInt * master_nodes_val = my_neighbor_list.master_nodes.values;
  
  for (UInt i = 0; i < nb_nodes_neigh; ++i) {
    test_output << " Impactor node: " << impact_nodes_val[i] << " has master nodes:";
    for(UInt mn = master_nodes_offset_val[i]; mn < master_nodes_offset_val[i+1]; ++mn) {
      test_output << " " << master_nodes_val[mn];
    }
    test_output << std::endl;
  }  

  delete my_contact;
  
  finalize();

  return EXIT_SUCCESS;
}
