/**
 * @file   test_single_spring_friction.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Thu May 26 17:03:01 2011
 *
 * @brief  test with analytical solution for friction of a single spring
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


#include <stdio.h>
#include <stdlib.h>

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
#include "unique_constant_fric_coef.hh"
#include "velocity_weakening_coulomb.hh"
#include "contact_neighbor_structure.hh"
#include "regular_grid_neighbor_structure.hh"
#include "contact_search.hh"
#include "contact_search_explicit.hh"

#ifdef AKANTU_USE_IOHELPER
#include <io_helper.hh>
#endif //AKANTU_USE_IOHELPER


using namespace akantu;

/* -------------------------------------------------------------------------- */

Surface impactor = 0;
Surface master = 1;
UInt the_node = 0;

const ElementType element_type = _triangle_3; 
#ifdef AKANTU_USE_IOHELPER
const iohelper::ElemType paraview_type = iohelper::TRIANGLE1; 
#endif //AKANTU_USE_IOHELPER
const char* mesh_name = "single_triangle.msh";
const char* folder_name = "single_spring_friction";
std::string result_name;
const bool viscous_mat = false;
const bool vel_weak_coulomb = false;
Real mu_s = 0.3;
Real mu_d = 0.2;
const UInt paraview_dump_int = 1e3;
const UInt file_dump_int = 1e2;
const UInt max_steps = 1e6;
const UInt stick_stop = 3;

Real n_k = 0.15;
Real normal_force;
Real start_displacement = 0.2;

Mesh * mesh;
SolidMechanicsModel * model;
UInt spatial_dimension = 2;
Real time_step_security = 3.;
Real time_step;
UInt nb_nodes;
UInt nb_elements;
UInt stick_counter = 0;

bool patch_test = false;
const Real tolerance = std::numeric_limits<Real>::epsilon() * 10.;

// variables of model
Real * coordinates; 
Real * displacement;
Real * current_position;
Real * velocity;
Real * residual;
Real * acceleration;
bool * boundary;
Real * force;
Real * mass;

std::map < std::string, ArrayBase* > restart_map;

#ifdef AKANTU_USE_IOHELPER
static void paraviewInit(iohelper::Dumper & dumper);
#endif //AKANTU_USE_IOHELPER

static void loadRestartInformation(ContactRigid * contact);
static void printPredictor(UInt step, 
			   std::ofstream & out_stream);
static void printCorrector(UInt step, 
			   ContactRigid * contact, 
			   std::ofstream & out_stream);
static void getStickInfo(ContactRigid * contact);
static bool testFloat(Real a, Real b, Real adm_error);

/* -------------------------------------------------------------------------- */
Int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  debug::setDebugLevel(dblWarning);

  // 1: name
  // 2: initial displacement
  // 3: friction coefficient
  // 4: ratio N/K with N normal force & K stiffness
  if (argc == 5) {
    result_name = argv[1];
    start_displacement = atof(argv[2]);
    mu_s = atof(argv[3]);
    n_k = atof(argv[4]);
  }
  else {
    result_name = "patch_test_output";
    patch_test = true;
  //    return EXIT_FAILURE;    
  }

  /// load mesh
  mesh = new Mesh(spatial_dimension);
  MeshIOMSH mesh_io;
  mesh_io.read(mesh_name, *mesh);

  /// build facet connectivity and surface id
  MeshUtils::buildFacets(*mesh);
  MeshUtils::buildSurfaceID(*mesh);

  /// declaration of model
  model = new SolidMechanicsModel(*mesh);

  nb_nodes = model->getFEM().getMesh().getNbNodes();
  nb_elements = model->getFEM().getMesh().getNbElement(element_type);
  std::cout << "Nb nodes : " << nb_nodes << " - nb elements : " << nb_elements << std::endl;

  /// model initialization
  model->initArrays();
  // initialize the vectors
  model->getForce().clear();
  model->getVelocity().clear();
  model->getAcceleration().clear();
  model->getDisplacement().clear();
  
  model->initExplicit();
  model->initModel();

  /// read and initialize material
  if(viscous_mat)
    model->readMaterials("material_elastic_caughey.dat");
  else
    model->readMaterials("material.dat");
  model->initMaterials();

  Real stable_time_step = model->getStableTimeStep();
  time_step = stable_time_step/time_step_security;
  model->setTimeStep(time_step);
  std::cout << "The time step is " << time_step << std::endl;

  // accessors to model elements
  coordinates      = mesh->getNodes().values;
  displacement     = model->getDisplacement().values;
  velocity         = model->getVelocity().values;
  acceleration     = model->getAcceleration().values;
  boundary         = model->getBoundary().values;
  residual         = model->getResidual().values;
  model->updateCurrentPosition();
  current_position = model->getCurrentPosition().values; 
  force            = model->getForce().values;
  mass             = model->getMass().values;
    
  model->assembleMassLumped();

  UInt nb_surfaces = mesh->getNbSurfaces();
  nb_surfaces += 1;

  /// contact declaration
  Contact * contact_structure = Contact::newContact(*model, 
						    _ct_rigid, 
						    _cst_expli, 
						    _cnst_regular_grid);
  ContactRigid * contact = dynamic_cast<ContactRigid *>(contact_structure);
  //contact->initContact(false);

  contact->addMasterSurface(master);
  contact->addImpactorSurfaceToMasterSurface(impactor, master);  

  /// define the friction law
  FrictionCoefficient *fric_coef;
  if(vel_weak_coulomb)
   fric_coef = new VelocityWeakeningCoulomb(*contact, master, mu_s, mu_d);
  else
    fric_coef = new UniqueConstantFricCoef(*contact, master, mu_s);

  // load restart file
  loadRestartInformation(contact);

  // set boundary condition
  displacement[the_node*spatial_dimension] = start_displacement;
  model->updateCurrentPosition();
  model->updateResidual();
  Real stiffness = std::abs(residual[the_node * spatial_dimension] / start_displacement);
  normal_force = n_k * stiffness;
  force[the_node*spatial_dimension+1] = -normal_force;

  if (start_displacement > tolerance) {
    std::cout << "Start displacement = " << start_displacement << " ; Residual = " << residual[the_node*spatial_dimension] << " -> Stiffness K = " << stiffness << std::endl;
    std::cout << "Mass = " << mass[the_node] << std::endl;
  }
  else
    std::cout << "No start displacement!! " << std::endl;

  /// initialize the paraview output
#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperParaview dumper;
#endif //AKANTU_USE_IOHELPER
  std::ofstream out_info;  
  if (!patch_test) {
    model->updateResidual();
#ifdef AKANTU_USE_IOHELPER
    paraviewInit(dumper);
#endif //AKANTU_USE_IOHELPER 
    
    /// output files
    std::stringstream name_info;
    name_info << "output_files/" << folder_name << "/global_info.dat";
    out_info.open(name_info.str().c_str());
    out_info << "%id time fnorm fres ffric ftot mu disp vel stick ekin epot etot" << std::endl; 
  }

  UInt step = 0;
  Real previous_vel = 0.;
  UInt count_cycle = 0;

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  while(stick_counter <= stick_stop && step <= max_steps) {
    
    // increase step
    step += 1;

    if(step % 1000 == 0) {
      std::cout << "passing step " << step << "/" << max_steps << "\r";
      std::cout.flush();
    }
    
    model->explicitPred();
    model->updateResidual();

    if (step % file_dump_int == 0 && !patch_test)
      printPredictor(step, out_info);

    contact->frictionPredictor();
    model->updateAcceleration();
    model->explicitCorr();
    contact->frictionCorrector();

    // see if node sticks
    getStickInfo(contact);
    
    // count numbers of cycles
    if (Math::are_float_equal(velocity[the_node * spatial_dimension], 0.) && 
   	!Math::are_float_equal(previous_vel, 0.))
      count_cycle++;
    previous_vel = velocity[the_node * spatial_dimension];

    if(step % file_dump_int == 0 && !patch_test)
      printCorrector(step, contact, out_info);

#ifdef AKANTU_USE_IOHELPER
    if(step % paraview_dump_int == 0 && !patch_test) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  }

  std::cout << "passing step " << step << "/" << max_steps << std::endl;

  Real final_displacement = displacement[the_node * spatial_dimension];

  // check patch test
  if(patch_test) {
    Real correct_final_disp = 0.0200181;
    Real admissible_error = 1e-07;    
    UInt correct_nb_cycle = 2;
    //if (!(Math::are_float_equal(final_displacement, correct_final_disp)) | !(count_cycle == correct_nb_cycle)) {
    if (!(testFloat(final_displacement, correct_final_disp, admissible_error)) | !(count_cycle == correct_nb_cycle)) {
      std::cout << "Final displacement is " << final_displacement << ", which should be " << correct_final_disp << std::endl;
      std::cout << "Final number of cycles is " << count_cycle << ", which should be " << correct_nb_cycle << std::endl;
      return EXIT_FAILURE;
    }
    else {
      std::cout << "Patch test successful!" << std::endl;
      return EXIT_SUCCESS;
    }
  }
  else {
    /// output files
    std::stringstream result_name_path;
    result_name_path << "output_files/" << folder_name << "/" << result_name << ".dat";
    std::ofstream out_result;
    out_result.open(result_name_path.str().c_str(), std::ios::app);
    if (!out_result.good()) {
      std::cout << "Could not open result file " << std::endl;
      return EXIT_FAILURE;
    }
    out_result << start_displacement << " " << final_displacement << " " << count_cycle << std::endl;
    out_result.close();
  }
  
  out_info.close();

  delete fric_coef;
  delete contact;
  delete model;
  delete mesh;
  
  finalize();

  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
void paraviewInit(iohelper::Dumper & dumper) {
  std::stringstream name;
  name << "paraview/" << folder_name << "/";

  dumper.SetMode(iohelper::TEXT);
  dumper.SetPoints(model->getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, "coordinates");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(element_type).values,
			 paraview_type, nb_elements, iohelper::C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().values,
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model->getVelocity().values,
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model->getResidual().values,
			  spatial_dimension, "force");
  dumper.AddNodeDataField(model->getMass().values,
			  spatial_dimension, "mass");
  dumper.AddNodeDataField(model->getForce().values,
			  spatial_dimension, "applied_force");
    
  dumper.AddElemDataField(model->getMaterial(0).getStrain(element_type).values, 
			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(model->getMaterial(0).getStress(element_type).values, 
			  spatial_dimension*spatial_dimension, "stress");
    
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("applied_force", 1);
  dumper.SetPrefix(name.str().c_str());
  dumper.Init();
  dumper.Dump();
}
#endif //AKANTU_USE_IOHELPER

/* -------------------------------------------------------------------------- */
void loadRestartInformation(ContactRigid * contact) {
    
  // boundary conditions
  Array<bool> * boundary_r = new Array<bool>(nb_nodes, spatial_dimension, false);
  // sticked nodes
  (*boundary_r)(1,0) = true;  (*boundary_r)(1,1) = true;
  (*boundary_r)(2,0) = true;  (*boundary_r)(2,1) = true;
  // the impactor node
  (*boundary_r)(0,1) = true;
  memcpy(boundary, boundary_r->values, spatial_dimension*nb_nodes*sizeof(bool));

  // set the active impactor node
  Array<bool> * ai_nodes = new Array<bool>(nb_nodes, 1, false);
  (*ai_nodes)(the_node) = true;
  restart_map["active_impactor_nodes"] = ai_nodes;

  // not defined master type, because won't use solve contact
  Array<ElementType> * et_nodes = new Array<ElementType>(nb_nodes, 1, _not_defined);
  restart_map["master_element_type"] = et_nodes;
  
  // master normal x=0; y=1
  Array<Real> * mn_nodes = new Array<Real>(0, spatial_dimension);
  Real normal[spatial_dimension];
  normal[0] = 0.; normal[1] = 1.;
  mn_nodes->push_back(normal);
  restart_map["master_normals"] = mn_nodes; 
  
  // node is sticking
  Array<bool> * is_nodes = new Array<bool>(nb_nodes, 2, false);
  (*is_nodes)(0,0) = true;  
  (*is_nodes)(0,1) = true;  
  restart_map["node_is_sticking"] = is_nodes;

  // no friction force
  Array<Real> * ff_nodes = new Array<Real>(0, spatial_dimension);
  Real force[spatial_dimension];
  force[0] = 0.; force[1] = 0.;
  ff_nodes->push_back(force);
  restart_map["friction_forces"] = ff_nodes;

  // the original stick position (for regularized friction)
  Array<Real> * sp_nodes = new Array<Real>(0, spatial_dimension);
  Real position[spatial_dimension];
  position[0] = 0.;
  position[1] = 0.;
  sp_nodes->push_back(position);
  restart_map["stick_positions"] = sp_nodes;

  // no residual forces
  Array<Real> * rf_nodes = new Array<Real>(0, spatial_dimension);
  Real r_force[spatial_dimension];
  r_force[0] = 0.; r_force[1] = 0.;
  rf_nodes->push_back(r_force);
  restart_map["residual_forces"] = rf_nodes;

  // no previous velocities
  Array<Real> * pv_nodes = new Array<Real>(0, spatial_dimension);
  Real vel[spatial_dimension];
  vel[0] = 0.; vel[1] = 0.;
  pv_nodes->push_back(vel);
  restart_map["previous_velocities"] = pv_nodes;

  // no friction strength
  Array<Real> * fr_nodes = new Array<Real>(0, 1);
  fr_nodes->push_back(0.);
  restart_map["friction_resistances"] = fr_nodes;

  contact->setRestartInformation(restart_map, master);
  delete boundary_r;
  delete ai_nodes;
  delete et_nodes;
  delete mn_nodes;
  delete is_nodes;
  delete ff_nodes;
  delete sp_nodes;
  delete rf_nodes;
  delete pv_nodes;
  delete fr_nodes;
}

/* -------------------------------------------------------------------------- */
void printPredictor(UInt step, std::ofstream & out_stream) {

  out_stream << step << " " << step*time_step << " " << residual[the_node*spatial_dimension + 1] << " " << residual[the_node * spatial_dimension];


}

/* -------------------------------------------------------------------------- */
void printCorrector(UInt step, ContactRigid * contact, std::ofstream & out_stream) {

  Real epot = model->getPotentialEnergy();
  Real ekin = model->getKineticEnergy();
  
  // find index of master surface in impactors_information 
  ContactRigid::SurfaceToImpactInfoMap::const_iterator it_imp;
  it_imp = contact->getImpactorsInformation().find(master);
  
  // find the total contact force and contact area
  ContactRigid::ImpactorInformationPerMaster * imp_info = it_imp->second;

  Real * friction_forces_val = imp_info->friction_forces->values;  
  bool * sticking_nodes_val = imp_info->node_is_sticking->values;
  
  out_stream << " " << friction_forces_val[0] << " " << residual[the_node * spatial_dimension] << " " << friction_forces_val[0] / residual[the_node * spatial_dimension + 1] << " " << displacement[the_node * spatial_dimension] << " " << velocity[the_node * spatial_dimension] << " " << sticking_nodes_val[0] << " " << ekin << " " << epot << " " << ekin+epot << std::endl;

  if (sticking_nodes_val[0])
    stick_counter++;
  else 
    stick_counter = 0;
}

/* -------------------------------------------------------------------------- */
void getStickInfo(ContactRigid * contact) {

  // find index of master surface in impactors_information 
  ContactRigid::SurfaceToImpactInfoMap::const_iterator it_imp;
  it_imp = contact->getImpactorsInformation().find(master);
  
  // find the total contact force and contact area
  ContactRigid::ImpactorInformationPerMaster * imp_info = it_imp->second;

  bool * sticking_nodes_val = imp_info->node_is_sticking->values;

  if (sticking_nodes_val[the_node*2])
    stick_counter++;
  else 
    stick_counter = 0;
}

/* -------------------------------------------------------------------------- */
bool testFloat(Real a, Real b, Real adm_error) {
  
  if (fabs(a-b) < adm_error)
    return true;
  else 
    return false;
}
