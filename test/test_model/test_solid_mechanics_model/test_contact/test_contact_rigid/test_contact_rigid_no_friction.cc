/**
 * @file   test_contact_rigid_no_friction.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Thu Oct 06 17:24:36 2011
 *
 * @brief  test rigid contact for all types
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
#include "contact_search.hh"
#include "contact_search_explicit.hh"

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper_tools.hh"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);

  const ElementType element_type = TYPE;
  UInt dim = Mesh::getSpatialDimension(element_type);

  /// load mesh
  Mesh my_mesh(dim);
  MeshIOMSH mesh_io;
  std::stringstream meshname_sstr; 
  meshname_sstr << element_type << ".msh";
  mesh_io.read(meshname_sstr.str().c_str(), my_mesh);

  /// build facet connectivity and surface id
  MeshUtils::buildFacets(my_mesh);
  MeshUtils::buildSurfaceID(my_mesh);

  UInt max_steps = 3; 
  unsigned int nb_nodes = my_mesh.getNbNodes();

  /// declaration of model
  SolidMechanicsModel  my_model(my_mesh);
  /// model initialization
  my_model.initVectors();
  // initialize the vectors
  my_model.getForce().clear();
  my_model.getVelocity().clear();
  my_model.getAcceleration().clear();
  my_model.getDisplacement().clear();

  Real * displacement = my_model.getDisplacement().values;
  Real * residual = my_model.getResidual().values;

  my_model.initExplicit();
  my_model.initModel();  
  my_model.readMaterials("material.dat");
  my_model.initMaterials();

  /// dump facet and surface information to paraview
#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperParaview dumper;
  paraviewInit(dumper, my_model, element_type, "para");
#endif //AKANTU_USE_IOHELPER

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

  /// define output file for testing
  std::stringstream filename_sstr; 
  filename_sstr << "test_contact_rigid_no_friction_" << element_type << ".out";
  std::ofstream test_output;
  test_output.open(filename_sstr.str().c_str());

  /*
  const NodesNeighborList & my_neighbor_list = dynamic_cast<const NodesNeighborList &>(my_contact->getContactSearch().getContactNeighborStructure(master).getNeighborList());

  UInt nb_nodes_neigh = my_neighbor_list.impactor_nodes.getSize();
  Vector<UInt> impact_nodes = my_neighbor_list.impactor_nodes;
  UInt * impact_nodes_val = impact_nodes.values;
  UInt * master_nodes_offset_val = my_neighbor_list.master_nodes_offset.values;
  UInt * master_nodes_val = my_neighbor_list.master_nodes.values;

  /// print impactor nodes
  test_output << "we have " << nb_nodes_neigh << " impactor nodes:";
  for (UInt i = 0; i < nb_nodes_neigh; ++i) {
    test_output << " " << impact_nodes_val[i];
  }
  test_output << std::endl;
  
  for (UInt i = 0; i < nb_nodes_neigh; ++i) {
    test_output << " Impactor node: " << impact_nodes_val[i] << " has master nodes:";
    for(UInt mn = master_nodes_offset_val[i]; mn < master_nodes_offset_val[i+1]; ++mn) {
      test_output << " " << master_nodes_val[mn];
    }
    test_output << std::endl;
  }
  */

  my_contact->initSearch(); // does nothing so far
  
  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(UInt s = 1; s <= max_steps; ++s) {

    test_output << std::endl << "passing step " << s << "/" << max_steps << std::endl;

    /// apply a displacement to the slave body
    if(s == 2) {
      Real * coord = my_mesh.getNodes().values;
      for(UInt n = 0; n < nb_nodes; ++n) {
	if(coord[n*dim + 0] > 0.5) {
	  displacement[n*dim+0] = -0.02;
	}
      }
    }

    /// integration
    my_model.explicitPred();
    my_model.initializeUpdateResidualData();

    /*
    /// compute the penetration list
    test_output << "Before solveContact" << std::endl;
    PenetrationList * my_penetration_list = new PenetrationList("penetration_list_1");
    const_cast<ContactSearch &>(my_contact->getContactSearch()).findPenetration(master, *my_penetration_list);

    UInt nb_nodes_pen = my_penetration_list->penetrating_nodes.getSize();
    Vector<UInt> pen_nodes = my_penetration_list->penetrating_nodes;
    UInt * pen_nodes_val = pen_nodes.values;
    test_output << "we have " << nb_nodes_pen << " penetrating nodes:" << std::endl;
    for (UInt i = 0; i < nb_nodes_pen; ++i) {
      test_output << "node " << pen_nodes_val[i] << " with disp:";
      for (UInt j=0; j<dim; ++j)
	test_output << " " << std::setprecision(10) << displacement[pen_nodes_val[i]*dim+j];
      test_output << std::endl;
    }
    test_output << std::endl;
    delete my_penetration_list;
    */

    /// solve contact
    my_contact->solveContact();

    /*
    /// compute the penetration list
    test_output << "After solveContact" << std::endl;
    PenetrationList * my_penetration_list_2 = new PenetrationList("penetration_list_2");
    const_cast<ContactSearch &>(my_contact->getContactSearch()).findPenetration(master, *my_penetration_list_2);

    UInt nb_nodes_pen_2 = my_penetration_list_2->penetrating_nodes.getSize();
    Vector<UInt> pen_nodes_2 = my_penetration_list_2->penetrating_nodes;
    UInt * pen_nodes_2_val = pen_nodes_2.values;
    test_output << "we have " << nb_nodes_pen_2 << " penetrating nodes:";
    for (UInt i = 0; i < nb_nodes_pen_2; ++i)
      test_output << " " << pen_nodes_2_val[i];
    test_output << std::endl;
    delete my_penetration_list_2;
    */

    /// compute the residual
    my_model.updateResidual(false);

    ContactRigid::SurfaceToImpactInfoMap::const_iterator it_imp;
    it_imp = my_contact->getImpactorsInformation().find(master);
    ContactRigid::ImpactorInformationPerMaster * imp_info = it_imp->second;
    UInt * active_imp_nodes_val = imp_info->active_impactor_nodes->values;
    test_output << "we have " << imp_info->active_impactor_nodes->getSize() << " active impactor nodes:" << std::endl;
    for (UInt i = 0; i < imp_info->active_impactor_nodes->getSize(); ++i) {
      UInt node = active_imp_nodes_val[i];
      test_output << "node " << node << " with disp:";
      for (UInt j=0; j<dim; ++j)
	test_output << " " << std::setprecision(10) << displacement[node*dim+j];
      test_output << " and force:";
      for (UInt j=0; j<dim; ++j)
	test_output << " " << std::setprecision(10) << residual[node*dim+j];
      test_output << std::endl;
    }

    // further integration
    my_contact->avoidAdhesion();
    my_model.updateAcceleration();
    my_model.explicitCorr();

#ifdef AKANTU_USE_IOHELPER
    paraviewDump(dumper);
#endif //AKANTU_USE_IOHELPER
  }

  delete my_contact;
  
  finalize();

  return EXIT_SUCCESS;
}
