/**
 * @file   test_contact_regular_grid.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Jan 17 11:13:42 2011
 *
 * @brief  test contact search for 2d case in explicit
 *
 * @section LICENSE
 *
 * <insert license here>
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
#include "contact_neighbor_structure.hh"
#include "regular_grid_neighbor_structure.hh"
#include "contact_search.hh"
#include "contact_search_3d_explicit.hh"

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  int dim = 2;
  const ElementType element_type = _triangle_3;

  /// load mesh
  Mesh my_mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("squares.msh", my_mesh);

  /// build facet connectivity and surface id
  MeshUtils::buildFacets(my_mesh,1,0);
  MeshUtils::buildSurfaceID(my_mesh);

  UInt max_steps = 3; 
  unsigned int nb_nodes = my_mesh.getNbNodes();

  /// dump facet and surface information to paraview
#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(TEXT);
  
  dumper.SetPoints(my_mesh.getNodes().values, dim, nb_nodes, "triangle_3_nodes_test-surface-extraction");
  dumper.SetConnectivity((int*)my_mesh.getConnectivity(_triangle_3).values,
   			 TETRA1, my_mesh.getNbElement(_triangle_3), C_MODE);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  /// declaration of model
  SolidMechanicsModel  my_model(my_mesh);
  /// model initialization
  my_model.initVectors();
  // initialize the vectors
  memset(my_model.getForce().values,        0, dim*nb_nodes*sizeof(Real));
  memset(my_model.getVelocity().values,     0, dim*nb_nodes*sizeof(Real));
  memset(my_model.getAcceleration().values, 0, dim*nb_nodes*sizeof(Real));
  memset(my_model.getDisplacement().values, 0, dim*nb_nodes*sizeof(Real));

  Real * displacement = my_model.getDisplacement().values;
  
  my_model.readMaterials("material.dat");
  my_model.initMaterials();
  my_model.initModel();

  Real time_step = my_model.getStableTimeStep();
  my_model.setTimeStep(time_step/10.);

  my_model.assembleMassLumped();

   /// contact declaration
  Contact * my_contact = Contact::newContact(my_model, 
					     _ct_rigid_no_fric, 
					     _cst_3d_expli, 
					     _cnst_regular_grid);

  my_contact->initContact(false);

  Surface master = 0;
  Surface impactor = 1;
  my_contact->addMasterSurface(master);
  
  my_model.updateCurrentPosition(); // neighbor structure uses current position for init
  my_contact->initNeighborStructure(master);

  const NodesNeighborList & my_neighbor_list = dynamic_cast<const NodesNeighborList &>(my_contact->getContactSearch().getContactNeighborStructure(master).getNeighborList());

  UInt nb_nodes_neigh = my_neighbor_list.impactor_nodes.getSize();
  Vector<UInt> impact_nodes = my_neighbor_list.impactor_nodes;
  UInt * impact_nodes_val = impact_nodes.values;

  /// print impactor nodes
  std::cout << "we have " << nb_nodes_neigh << " impactor nodes:";
  for (UInt i = 0; i < nb_nodes_neigh; ++i) {
    std::cout << " " << impact_nodes_val[i];
  }
  std::cout << std::endl;

  UInt * master_nodes_offset_val = my_neighbor_list.master_nodes_offset.values;
  UInt * master_nodes_val = my_neighbor_list.master_nodes.values;
  
  for (UInt i = 0; i < nb_nodes_neigh; ++i) {
    std::cout << " Impactor node: " << impact_nodes_val[i] << " has master nodes:";
    for(UInt mn = master_nodes_offset_val[i]; mn < master_nodes_offset_val[i+1]; ++mn) {
      std::cout << " " << master_nodes_val[mn];
    }
    std::cout << std::endl;
  }

  my_contact->initSearch(); // does nothing so far
  
  std::cout << std::endl << "epsilon = " << std::numeric_limits<Real>::epsilon() << std::endl;

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(UInt s = 1; s <= max_steps; ++s) {

    std::cout << std::endl << "passing step " << s << "/" << max_steps << std::endl;

    /// apply a displacement to the slave body
    if(s == 2) {
      Real * coord = my_mesh.getNodes().values;
      for(UInt n = 0; n < nb_nodes; ++n) {
	if(coord[n*dim + 0] > 1.0) {
	  displacement[n*dim+0] = -0.02;
	}
      }
      /*
      UInt nb_elements = my_mesh.getNbElement(element_type);
      UInt nb_nodes_element = my_mesh.getNbNodesPerElement(element_type);
      Vector<UInt> element_mat = my_model.getElementMaterial(element_type);
      UInt * element_mat_val = element_mat.values;
      UInt * connectivity = my_mesh.getConnectivity(element_type).values;
      for(UInt el = 0; el < nb_elements; ++el) {
	std::cout << "element: " << el << " with mat: " <<  element_mat_val[el] << std::endl;
	if(element_mat_val[el] == impactor) {
	  for(UInt n = 0; n < nb_nodes_element; ++n) {
	    displacement[connectivity[el * nb_nodes_element + n]+2] = -0.2;
	  }
	}
	}*/
    }

    /// central difference predictor
    my_model.explicitPred();
    
    /// update current position
    my_model.updateCurrentPosition();

    /// compute the penetration list
    std::cout << "Before solveContact" << std::endl;
    PenetrationList * my_penetration_list = new PenetrationList();
    const_cast<ContactSearch &>(my_contact->getContactSearch()).findPenetration(master, *my_penetration_list);

    UInt nb_nodes_pen = my_penetration_list->penetrating_nodes.getSize();
    Vector<UInt> pen_nodes = my_penetration_list->penetrating_nodes;
    UInt * pen_nodes_val = pen_nodes.values;
    std::cout << "we have " << nb_nodes_pen << " penetrating nodes:";
    for (UInt i = 0; i < nb_nodes_pen; ++i)
      std::cout << " " << pen_nodes_val[i];
    std::cout << std::endl;
    delete my_penetration_list;

    /// solve contact
    my_contact->solveContact();

    /// compute the penetration list
    std::cout << "After solveContact" << std::endl;
    PenetrationList * my_penetration_list_2 = new PenetrationList();
    const_cast<ContactSearch &>(my_contact->getContactSearch()).findPenetration(master, *my_penetration_list_2);

    UInt nb_nodes_pen_2 = my_penetration_list_2->penetrating_nodes.getSize();
    Vector<UInt> pen_nodes_2 = my_penetration_list_2->penetrating_nodes;
    UInt * pen_nodes_2_val = pen_nodes_2.values;
    std::cout << "we have " << nb_nodes_pen_2 << " penetrating nodes:";
    for (UInt i = 0; i < nb_nodes_pen_2; ++i)
      std::cout << " " << pen_nodes_2_val[i];
    std::cout << std::endl;
    delete my_penetration_list_2;

    /// compute the residual
    my_model.updateResidual(false);
    
    /// compute the acceleration
    my_model.updateAcceleration();

    /// central difference corrector
    my_model.explicitCorr();
  }

  delete my_contact;
  
  finalize();

  return EXIT_SUCCESS;
}
