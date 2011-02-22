/**
 * @file   test_contact_rigid_ball_2d.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Feb 22 09:04:42 2011
 *
 * @brief  test rigid contact for a 2d ball
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
  const UInt paraview_type = TRIANGLE1;

  UInt max_steps = 200000;

  /// load mesh
  Mesh my_mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("ball_2d.msh", my_mesh);

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

  Surface impactor = 4;
  Surface rigid_body_surface[4];
  rigid_body_surface[0] = 0;
  rigid_body_surface[1] = 1;
  rigid_body_surface[2] = 2;
  rigid_body_surface[3] = 3;

  Surface master = 5;
  //master[0] = 5;
  //master[1] = 6;
  //master[2] = 7;
  //master[3] = 8;

  // modify surface id
  UInt nb_surfaces = my_mesh.getNbSurfaces();
  my_mesh.setNbSurfaces(++nb_surfaces); 
  ElementType surface_element_type = my_mesh.getFacetElementType(element_type);
  UInt nb_surface_element = my_model.getFEM().getMesh().getNbElement(surface_element_type);
  UInt * surface_id_val = my_mesh.getSurfaceId(surface_element_type).values;
  for(UInt i=0; i < nb_surface_element; ++i) {
    if ((surface_id_val[i] == rigid_body_surface[0]) ||
	(surface_id_val[i] == rigid_body_surface[1]) ||
	(surface_id_val[i] == rigid_body_surface[2]) ||
	(surface_id_val[i] == rigid_body_surface[3])) {
      Real barycenter[dim];
      Real * barycenter_p = &barycenter[0];
      my_mesh.getBarycenter(i,surface_element_type,barycenter_p);
      if((barycenter_p[0] > -0.001) && 
	 (barycenter_p[0] <  0.501) &&
	 (barycenter_p[1] > -0.001) &&
	 (barycenter_p[1] <  0.701)) {
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

  my_contact->addMasterSurface(master);
  my_contact->addImpactorSurfaceToMasterSurface(impactor, master);

  my_model.updateCurrentPosition(); // neighbor structure uses current position for init
  my_contact->initNeighborStructure(master);
  my_contact->initSearch(); // does nothing so far

  // boundary conditions
  Real * coordinates = my_mesh.getNodes().values;
  bool * boundary = my_model.getBoundary().values;
  Real * velocity_val = my_model.getVelocity().values;
  for (UInt i=0; i < nb_nodes; ++i) {
    Real x_coord = coordinates[i*dim];
    Real y_coord = coordinates[i*dim + 1];
    if((x_coord > 0) &&
       (x_coord < 0.5) &&
       (y_coord > 0) &&
       (y_coord < 0.7)) {
      velocity_val[i*dim]   = 100;
      velocity_val[i*dim+1] = -100;
    }
  }  

  
#ifdef AKANTU_USE_IOHELPER
  /// initialize the paraview output
  DumperParaview dumper;
  dumper.SetMode(TEXT);
  dumper.SetPoints(my_model.getFEM().getMesh().getNodes().values, dim, nb_nodes, "coordinates_2d");
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
  dumper.SetPrefix("paraview/ball_2d/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  std::ofstream energy;
  energy.open("ball_2d_energy.csv");
  energy << "%id,kin,pot,tot" << std::endl;


  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(UInt s = 1; s <= max_steps; ++s) {

    if(s % 10 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;

    if(s % 2000 == 0){
      my_model.updateCurrentPosition();
      my_contact->updateContact();    
    }

    my_model.explicitPred();
   
    my_model.initializeUpdateResidualData();

    my_contact->solveContact();

    my_model.updateResidual(false);

    my_contact->avoidAdhesion();

    my_contact->addFriction();

    my_model.updateAcceleration();
    my_model.explicitCorr();

    my_contact->addSticking();

    Real epot = my_model.getPotentialEnergy();
    Real ekin = my_model.getKineticEnergy();
    energy << s << "," << ekin << "," << epot << "," << ekin+epot << std::endl;


#ifdef AKANTU_USE_IOHELPER
    if(s % 100 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  }

  energy.close();

  delete my_contact;
 
  finalize();

  return EXIT_SUCCESS;
}
