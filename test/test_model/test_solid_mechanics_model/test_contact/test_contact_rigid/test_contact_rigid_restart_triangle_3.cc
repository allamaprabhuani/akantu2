/**
 * @file   test_contact_rigid_restart_triangle_3.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Fri Apr 29 11:19:46 2011
 *
 * @brief  test restart functions for contact rigid
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
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  UInt dim = 2;
  const ElementType element_type = _triangle_3;

  akantu::initialize(&argc, &argv);

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
  
  dumper.SetPoints(my_mesh.getNodes().values, dim, nb_nodes, "triangle_3_restart_test");
  dumper.SetConnectivity((int*)my_mesh.getConnectivity(element_type).values,
   			 TRIANGLE1, my_mesh.getNbElement(element_type), C_MODE);
  dumper.SetPrefix("paraview/triangle/");
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

  my_model.initExplicit();
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
    my_model.initializeUpdateResidualData();

    /// compute the penetration list
    std::cout << "Before solveContact" << std::endl;
    PenetrationList * my_penetration_list = new PenetrationList("penetration_list_1");
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
    PenetrationList * my_penetration_list_2 = new PenetrationList("penetration_list_2");
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

    my_contact->avoidAdhesion();
    
    /// compute the acceleration
    my_model.updateAcceleration();

    /// central difference corrector
    my_model.explicitCorr();
  }

  std::map < std::string, VectorBase* > restart_map;
  my_contact->getRestartInformation(restart_map, master);
  
  std::map < std::string, VectorBase* >::iterator it;
  
  it = restart_map.find("active_impactor_nodes");
  Vector<bool> * ai_nodes = (Vector<bool> *)(it->second);
  if(it == restart_map.end()) {
    std::cout << "could not find map entry for active impactor nodes" << std::endl;
  }

  it = restart_map.find("master_element_type");
  Vector<ElementType> * et_nodes = (Vector<ElementType> *)(it->second);
  if(it == restart_map.end()) {
    std::cout << "could not find map entry master element type" << std::endl;
  }

  it = restart_map.find("master_normals");
  Vector<Real> * mn_nodes = (Vector<Real> *)(it->second);
  if(it == restart_map.end()) {
    std::cout << "could not find map entry for master normals" << std::endl;
  }

  it = restart_map.find("node_is_sticking");
  Vector<bool> * is_nodes = (Vector<bool> *)(it->second);
  if(it == restart_map.end()) {
    std::cout << "could not find map entry node is sticking" << std::endl;
  }

  it = restart_map.find("friction_forces");
  Vector<Real> * ff_nodes = (Vector<Real> *)(it->second);
  if(it == restart_map.end()) {
    std::cout << "could not find map entry friction forces" << std::endl;
  }


  std::cout << "Active impactor nodes in the restart map:" << std::endl;
  for (UInt i=0; i<nb_nodes; ++i) {
    if ((*ai_nodes)(i)) {
      std::cout << "node: " << i << ", master element type: " << (*et_nodes)(i) << std::endl;
      for (UInt d=0; d<dim; ++d) {
	std::cout << "master normal direction = " << d << ": " << (*mn_nodes)(i,d);
	std::cout << "friction force direction = " << d << ": " << (*ff_nodes)(i,d) << std::endl;
      }
      for (UInt j=0; j<2; ++j) {
	std::cout << "stick " << j << ": " << (*is_nodes)(i,j);
      }
      std::cout << std::endl;
    }
  }

  Contact * contact_restart = Contact::newContact(my_model, 
						  _ct_rigid, 
						  _cst_expli, 
						  _cnst_regular_grid,
						  "contact_restart");

  ContactRigid * my_contact_restart = dynamic_cast<ContactRigid *>(contact_restart);
  
  my_contact_restart->initContact(false);
  my_contact_restart->addMasterSurface(master);
  my_contact_restart->addImpactorSurfaceToMasterSurface(impactor, master);
  my_contact_restart->setRestartInformation(restart_map, master);

  ContactRigid::SurfaceToImpactInfoMap::const_iterator it_imp;
  it_imp = my_contact_restart->getImpactorsInformation().find(master);
  ContactRigid::ImpactorInformationPerMaster * impactor_info = it_imp->second;

  UInt * active_nodes = impactor_info->active_impactor_nodes->values;
  ElementType * element_type_imp = &(*impactor_info->master_element_type)[0];  
  Real * master_normal = impactor_info->master_normals->values;
  bool * node_stick = impactor_info->node_is_sticking->values;
  Real * friction_force = impactor_info->friction_forces->values;

  UInt nb_active_nodes = impactor_info->active_impactor_nodes->getSize();
  std::cout << std::endl << "Active impactor nodes in restart contact:" << std::endl;
  for (UInt i=0; i<nb_active_nodes; ++i, ++active_nodes, ++element_type_imp) {
    std::cout << "node: " <<  (*active_nodes) << ", master element type: " << *element_type_imp << std::endl;
    for (UInt d=0; d<dim; ++d, ++master_normal, ++friction_force) {
	std::cout << "master normal direction = " << d << ": " << *master_normal;
	std::cout << "friction force direction = " << d << ": " << *friction_force << std::endl;
      }
    for (UInt j=0; j<2; ++j, ++node_stick) {
      std::cout << "stick " << j << ": " << *node_stick;
    }
	   std::cout << std::endl;
  }

  delete my_contact;
  delete my_contact_restart;
  
  finalize();

  return EXIT_SUCCESS;
}
