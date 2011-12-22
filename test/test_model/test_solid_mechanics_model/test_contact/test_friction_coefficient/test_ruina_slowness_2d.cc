/**
 * @file   test_ruina_slowness_2d.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Mar 07 15:56:42 2011
 *
 * @brief  the slowness law of the simplified dieterich friction coefficient
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
#include "friction_coefficient.hh"
#include "simplified_dieterich_fric_coef.hh"
#include "ruina_slowness_fric_coef.hh"
#include "unique_constant_fric_coef.hh"
#include "contact_neighbor_structure.hh"
#include "regular_grid_neighbor_structure.hh"
#include "contact_search.hh"
#include "contact_search_explicit.hh"

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);

  UInt dim = 2;
  const ElementType element_type = _triangle_3;
#ifdef AKANTU_USE_IOHELPER
  const UInt paraview_type = iohelper::TRIANGLE1;
#endif //AKANTU_USE_IOHELPER
  
  //UInt max_steps = 200000;
  UInt imposing_steps = 1000;
  Real max_displacement = -0.01;

  UInt damping_steps = 100000;
  UInt damping_interval = 200;
  Real damping_ratio = 0.99;

  Real sliding_velocity = 10.;
  UInt nb_inc_vel_steps = 20000;
  Real updated_displacement = 0.;

  UInt additional_steps = 200000;

  UInt max_steps = imposing_steps + damping_steps + nb_inc_vel_steps + additional_steps;

  /// load mesh
  Mesh my_mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("sliding_cube_2d.msh", my_mesh);

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

  my_model.initExplicit();
  my_model.initModel();
  my_model.readMaterials("material.dat");
  my_model.initMaterials();

  UInt nb_element = my_model.getFEM().getMesh().getNbElement(element_type);

  Real time_step = my_model.getStableTimeStep();
  my_model.setTimeStep(time_step/10.);

  //UInt max_steps = static_cast<UInt>(needed_time*10/time_step) + imposing_steps + damping_steps + additional_steps;
  std::cout << "The number of time steps is found to be: " << max_steps << std::endl;
  Real * velocity_val = my_model.getVelocity().values;
  
  my_model.assembleMassLumped();

  Surface impactor = 0;
  Surface rigid_body_surface = 1;
  Surface master = 2;

  // modify surface id
  UInt nb_surfaces = my_mesh.getNbSurfaces();
  my_mesh.setNbSurfaces(++nb_surfaces); 
  ElementType surface_element_type = my_mesh.getFacetElementType(element_type);
  UInt nb_surface_element = my_model.getFEM().getMesh().getNbElement(surface_element_type);
  UInt * surface_id_val = my_mesh.getSurfaceID(surface_element_type).values;
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

  my_contact->addMasterSurface(master);
  my_contact->addImpactorSurfaceToMasterSurface(impactor, master);  

  /*
  UniqueConstantFricCoef * fric_coef = new UniqueConstantFricCoef(*my_contact, master);
  fric_coef->setParam("mu", "0.2");
*/

  RuinaSlownessFricCoef<true> * fric_coef = new RuinaSlownessFricCoef<true>(*my_contact, master);
  fric_coef->setParam("mu_zero", "0.2");
  fric_coef->setParam("a_factor", "0.002");
  fric_coef->setParam("b_factor", "0.08");
  fric_coef->setParam("v_normalizer", "0.0001");
  fric_coef->setParam("theta_normalizer", "0.002");
  fric_coef->setParam("d_zero", "0.00001");
  //my_contact->setFrictionCoefficient(fric_coef);
  

  my_model.updateCurrentPosition(); // neighbor structure uses current position for init
  my_contact->initNeighborStructure(master);
  my_contact->initSearch(); // does nothing so far

  // boundary conditions
  Vector<UInt> * top_nodes = new Vector<UInt>(0, 1);
  Vector<UInt> * push_nodes = new Vector<UInt>(0, 1);
  Real * coordinates = my_mesh.getNodes().values;
  Real * displacement = my_model.getDisplacement().values;
  bool * boundary = my_model.getBoundary().values;
  UInt * surface_to_nodes_offset = my_contact->getSurfaceToNodesOffset().values;
  UInt * surface_to_nodes        = my_contact->getSurfaceToNodes().values;
  // normal force boundary conditions
  for(UInt n = surface_to_nodes_offset[impactor]; n < surface_to_nodes_offset[impactor+1]; ++n) {
    UInt node = surface_to_nodes[n];
    Real x_coord = coordinates[node*dim];
    Real y_coord = coordinates[node*dim + 1];
    if (y_coord > -0.00001) {
      boundary[node*dim + 1] = true;
      top_nodes->push_back(node);
    }
    if (x_coord < 0.00001) {
      boundary[node*dim] = true;
      push_nodes->push_back(node);
    }
  }
  // ground boundary conditions
  for(UInt n = surface_to_nodes_offset[rigid_body_surface]; n < surface_to_nodes_offset[rigid_body_surface+1]; ++n) {
    UInt node = surface_to_nodes[n];
    Real y_coord = coordinates[node*dim + 1];
    if (y_coord < -1.49)
      boundary[node*dim]     = true;
      boundary[node*dim + 1] = true;
  }
  UInt * top_nodes_val = top_nodes->values;
  UInt * push_nodes_val = push_nodes->values;

  my_model.updateResidual();
#ifdef AKANTU_USE_IOHELPER
  /// initialize the paraview output
  iohelper::DumperParaview dumper;
  dumper.SetMode(iohelper::TEXT);
  dumper.SetPoints(my_model.getFEM().getMesh().getNodes().values, dim, nb_nodes, "coord_ruina_slowness_2d");
  dumper.SetConnectivity((int *)my_model.getFEM().getMesh().getConnectivity(element_type).values,
			 paraview_type, nb_element, iohelper::C_MODE);
  dumper.AddNodeDataField(my_model.getDisplacement().values,
			  dim, "displacements");
  dumper.AddNodeDataField(my_model.getVelocity().values, dim, "velocity");
  dumper.AddNodeDataField(my_model.getResidual().values, dim, "force");
  dumper.AddNodeDataField(my_model.getForce().values, dim, "applied_force");
  dumper.AddElemDataField(my_model.getMaterial(0).getStrain(element_type).values, dim*dim, "strain");
  dumper.AddElemDataField(my_model.getMaterial(0).getStress(element_type).values, dim*dim, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("applied_force", 1);
  dumper.SetPrefix("paraview/ruina_slowness_2d/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  std::ofstream out_info;
  out_info.open("ruina_slowness_2d.csv");
  out_info << "%id,ftop,fcont,zone,stickNode,contNode" << std::endl;

  std::ofstream energy;
  energy.open("ruina_slowness_2d_energy.csv");
  energy << "%id,kin,pot,tot" << std::endl;

  Real * current_position = my_model.getCurrentPosition().values; 

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(UInt s = 1; s <= max_steps; ++s) {

    if(s % 10 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;

    if(s % 20000 == 0){
      my_model.updateCurrentPosition();
      my_contact->updateContact();    
    }

    // impose normal displacement
    if(s <= imposing_steps) {
      Real current_displacement = max_displacement/(static_cast<Real>(imposing_steps))*s;
      for(UInt n=0; n<top_nodes->getSize(); ++n) {
	UInt node = top_nodes_val[n];
	displacement[node*dim + 1] = current_displacement;
      }
    }
    
    // damp velocity in order to find equilibrium
    if(s > imposing_steps && s < imposing_steps+damping_steps && s%damping_interval == 0) {
      for (UInt i=0; i < nb_nodes; ++i) {
	for (UInt j=0; j < dim; ++j)
	  velocity_val[i*dim + j] *= damping_ratio;
      }
    }

    // give initial velocity
    if(s > imposing_steps+damping_steps) {
      Real t_step = my_model.getTimeStep();
      Real additional_disp;
      if (s > imposing_steps+damping_steps+nb_inc_vel_steps)
	additional_disp = sliding_velocity * t_step;
      else
	additional_disp = (s - (imposing_steps+damping_steps))/nb_inc_vel_steps*sliding_velocity*t_step;
      updated_displacement += additional_disp;
      for(UInt n=0; n<push_nodes->getSize(); ++n) {
	UInt node = push_nodes_val[n];
	displacement[node*dim] = updated_displacement;
      }
    }

    my_model.explicitPred();
   
    my_model.initializeUpdateResidualData();

    my_contact->solveContact();

    my_model.updateResidual(false);

    my_contact->avoidAdhesion();

    my_contact->frictionPredictor();

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
    Real contact_force = 0.;
    Real contact_zone = 0.;
    for (UInt i = 0; i < imp_info->active_impactor_nodes->getSize(); ++i) {
      UInt node = active_imp_nodes_val[i];
      contact_force += residual[node*dim + 1];
      contact_zone = std::max(contact_zone, current_position[node*dim]);
    }

    out_info << s << "," << top_force << "," << contact_force << "," << contact_zone << ",";

    my_model.updateAcceleration();

    const Vector<bool> * sticking_nodes = imp_info->node_is_sticking;
    bool * sticking_nodes_val = sticking_nodes->values;
    UInt nb_sticking_nodes = 0;
    for (UInt i = 0; i < imp_info->active_impactor_nodes->getSize(); ++i) {
      if(sticking_nodes_val[i*2])
	nb_sticking_nodes++;
    }

    out_info << nb_sticking_nodes << "," << imp_info->active_impactor_nodes->getSize() << std::endl;

    my_model.explicitCorr();

    my_contact->frictionCorrector();

    Real epot = my_model.getPotentialEnergy();
    Real ekin = my_model.getKineticEnergy();
    energy << s << "," << ekin << "," << epot << "," << ekin+epot << std::endl;


#ifdef AKANTU_USE_IOHELPER
    if(s % 100 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  }

  out_info.close();
  energy.close();

  delete fric_coef;
  delete my_contact;
 
  finalize();

  return EXIT_SUCCESS;
}
