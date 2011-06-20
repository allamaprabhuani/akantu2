/**
 * @file   test_contact_rigid_no_friction_force_3d.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Jan 24 11:56:42 2011
 *
 * @brief  test force for rigid contact 3d in explicit
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
  UInt dim = 3;
  const ElementType element_type = _tetrahedron_4;
#ifdef AKANTU_USE_IOHELPER
  const UInt paraview_type = TETRA1;
#endif //AKANTU_USE_IOHELPER
  
  UInt imposing_steps = 1000;
  Real max_displacement = -0.01;

  UInt damping_steps = 80000;
  UInt damping_interval = 50;
  Real damping_ratio = 0.985;

  UInt max_steps = damping_steps;

  /// load mesh
  Mesh my_mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("force_3d.msh", my_mesh);

  /// build facet connectivity and surface id
  MeshUtils::buildFacets(my_mesh,1,0);
  MeshUtils::buildSurfaceID(my_mesh);

  UInt nb_nodes = my_mesh.getNbNodes();
  
  /// declaration of model
  SolidMechanicsModel  my_model(my_mesh);
  /// model initialization
  my_model.initVectors();
  // initialize the vectors
  memset(my_model.getForce().values,        0,     dim*nb_nodes*sizeof(Real));
  memset(my_model.getVelocity().values,     0,     dim*nb_nodes*sizeof(Real));
  memset(my_model.getAcceleration().values, 0,     dim*nb_nodes*sizeof(Real));
  memset(my_model.getDisplacement().values, 0,     dim*nb_nodes*sizeof(Real));
  memset(my_model.getBoundary().values,     false, dim*nb_nodes*sizeof(bool));

  my_model.initModel();
  my_model.readMaterials("material.dat");
  my_model.initMaterials();

  UInt nb_element = my_model.getFEM().getMesh().getNbElement(element_type);

  Real time_step = my_model.getStableTimeStep();
  my_model.setTimeStep(time_step/10.);

  my_model.assembleMassLumped();

  Real * velocity_val = my_model.getVelocity().values;

  Surface impactor = 0;
  Surface rigid_body_surface = 1;
  Surface master = 2;

  // modify surface id
  UInt nb_surfaces = my_mesh.getNbSurfaces();
  my_mesh.setNbSurfaces(++nb_surfaces); 
  ElementType surface_element_type = my_mesh.getFacetElementType(element_type);
  //UInt * connectivity = my_mesh.getConnectivity(surface_element_type).values;
  //UInt nb_nodes_elem = Mesh::getNbNodesPerElement(surface_element_type);
  UInt nb_surface_element = my_model.getFEM().getMesh().getNbElement(surface_element_type);
  UInt * surface_id_val = my_mesh.getSurfaceId(surface_element_type).values;
  for(UInt i=0; i < nb_surface_element; ++i) {
    if (surface_id_val[i] == rigid_body_surface) {
      Real barycenter[dim];
      Real * barycenter_p = &barycenter[0];
      my_mesh.getBarycenter(i,surface_element_type,barycenter_p);
      if(barycenter_p[1] > -1.001) {
	surface_id_val[i] = master;
      }
    }
  }

   /// contact declaration
  Contact * contact = Contact::newContact(my_model, 
					  _ct_rigid, 
					  _cst_expli, 
					  _cnst_regular_grid);

  ContactRigid * my_contact = dynamic_cast<ContactRigid *>(contact);

  my_contact->initContact(false);

  //  Surface master = 1;
  my_contact->addMasterSurface(master);
  my_contact->addImpactorSurfaceToMasterSurface(impactor, master);

  my_model.updateCurrentPosition(); // neighbor structure uses current position for init
  my_contact->initNeighborStructure(master);
  my_contact->initSearch(); // does nothing so far

  // boundary conditions
  Vector<UInt> * top_nodes = new Vector<UInt>(0, 1);
  Real * coordinates = my_mesh.getNodes().values;
  Real * displacement = my_model.getDisplacement().values;
  bool * boundary = my_model.getBoundary().values;
  UInt * surface_to_nodes_offset = my_contact->getSurfaceToNodesOffset().values;
  UInt * surface_to_nodes        = my_contact->getSurfaceToNodes().values;
  // symetry boundary conditions
  for(UInt n = surface_to_nodes_offset[impactor]; n < surface_to_nodes_offset[impactor+1]; ++n) {
    UInt node = surface_to_nodes[n];
    Real y_coord = coordinates[node*dim + 1];
    if (y_coord > -0.00001) {
      boundary[node*dim + 1] = true;
      top_nodes->push_back(node);
    }
  }
  // ground boundary conditions
  for(UInt n = surface_to_nodes_offset[rigid_body_surface]; n < surface_to_nodes_offset[rigid_body_surface+1]; ++n) {
    UInt node = surface_to_nodes[n];
    Real y_coord = coordinates[node*dim + 1];
    if (y_coord < -1.19999) {
      boundary[node*dim]     = true;
      boundary[node*dim + 1] = true;
      boundary[node*dim + 2] = true;
    }
  }
  UInt * top_nodes_val = top_nodes->values;
  
  my_model.updateResidual();
#ifdef AKANTU_USE_IOHELPER
  /// initialize the paraview output
  DumperParaview dumper;
  dumper.SetMode(TEXT);
  dumper.SetPoints(my_model.getFEM().getMesh().getNodes().values, dim, nb_nodes, "coordinates_force_3d");
  dumper.SetConnectivity((int *)my_model.getFEM().getMesh().getConnectivity(element_type).values,
			 paraview_type, nb_element, C_MODE);
  dumper.AddNodeDataField(my_model.getDisplacement().values,
			  dim, "displacements");
  dumper.AddNodeDataField(my_model.getVelocity().values, dim, "velocity");
  dumper.AddNodeDataField(my_model.getResidual().values, dim, "force");
  dumper.AddNodeDataField(my_model.getForce().values, dim, "applied_force");
  dumper.AddElemDataField(my_model.getMaterial(0).getStrain(element_type).values, dim*dim, "strain");
  dumper.AddElemDataField(my_model.getMaterial(0).getStress(element_type).values, dim*dim, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("applied_force", 1);
  dumper.SetPrefix("paraview/force_3d/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  std::ofstream force_out;
  force_out.open("output_files/force_3d.csv");
  force_out << "%id,ftop,fcont,zone" << std::endl;


  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(UInt s = 1; s <= max_steps; ++s) {

    if(s % 10 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;

    if(s == imposing_steps){
      my_model.updateCurrentPosition();
      my_contact->updateContact();    
    }
    
    if(s <= imposing_steps) {
      Real current_displacement = max_displacement/(static_cast<Real>(imposing_steps))*s;
      for(UInt n=0; n<top_nodes->getSize(); ++n) {
	UInt node = top_nodes_val[n];
	displacement[node*dim + 1] = current_displacement;
      }
    }

    // damp velocity in order to find equilibrium
    if(s < damping_steps && s%damping_interval == 0) {
      for (UInt i=0; i < nb_nodes; ++i) {
	for (UInt j=0; j < dim; ++j)
	  velocity_val[i*dim + j] *= damping_ratio;
      }
    }
    
    my_model.explicitPred();
   
    my_model.initializeUpdateResidualData();

    my_contact->solveContact();

    my_model.updateResidual(false);

    // find the total force applied at the imposed displacement surface (top) 
    Real * residual = my_model.getResidual().values; 
    Real top_force = 0.;
    for(UInt n=0; n<top_nodes->getSize(); ++n) {
      UInt node = top_nodes_val[n];
      top_force += residual[node*dim + 1];
    }

    // find index of master surface in impactors_information 
    ContactRigid::SurfaceToImpactInfoMap::const_iterator it_imp;
    it_imp = my_contact->getImpactorsInformation().find(master);

    // find the total contact force and contact area
    ContactRigid::ImpactorInformationPerMaster * imp_info = it_imp->second;
    UInt * active_imp_nodes_val = imp_info->active_impactor_nodes->values;
    Real * current_position = my_model.getCurrentPosition().values; 
    Real contact_force = 0.;
    Real contact_zone = 0.;
    for (UInt i = 0; i < imp_info->active_impactor_nodes->getSize(); ++i) {
      UInt node = active_imp_nodes_val[i];
      contact_force += residual[node*dim + 1];
      contact_zone = std::max(contact_zone, current_position[node*dim]); 
    }

    force_out << s << "," << top_force << "," << contact_force << "," << contact_zone << std::endl;

    // test if correct result is found
    if((s == max_steps) && (std::abs(top_force - 2.1e+09) > 1e4)) {
      std::cout << "top_force = " << top_force << " but should be = 2.1e+09" << std::endl;
      return EXIT_FAILURE;
    }

    my_model.updateAcceleration();
    my_model.explicitCorr();

#ifdef AKANTU_USE_IOHELPER
    if(s % 1000 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  }

  force_out.close();

  delete my_contact;
 
  finalize();

  return EXIT_SUCCESS;
}
