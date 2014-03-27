/**
 * @file   test_contact_rigid_restart.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Mon May 02 11:22:55 2011
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
#include "friction_coefficient.hh"
#include "unique_constant_fric_coef.hh"
#include "contact_neighbor_structure.hh"
#include "regular_grid_neighbor_structure.hh"
#include "contact_search.hh"
#include "contact_search_explicit.hh"

#include <io_helper.hh>
#include <reader_restart.hh>

using namespace akantu;

static void restartReaderInit(iohelper::ReaderRestart & reader);
static void loadRestartInformation(ContactRigid * contact, std::map < std::string, ArrayBase* > & map);
static void DumpRestart(SolidMechanicsModel  & my_model, std::map < std::string, ArrayBase* > & map);
static void printRestartMap(std::map < std::string, ArrayBase* > & map);
static void printContact(ContactRigid * contact);
static void freeMap(std::map < std::string, ArrayBase* > & map);

UInt dim = 2;
const ElementType element_type = _triangle_3;
const iohelper::ElemType paraview_type = iohelper::TRIANGLE1;

Surface master;
UInt nb_nodes;
UInt nb_elements;
UInt nb_components;
Real * displacement;
bool * boundary;

std::ofstream test_output;

int main(int argc, char *argv[])
{
  /// define output file for testing
  std::stringstream filename_sstr;
  filename_sstr << "test_contact_rigid_restart.out";
  test_output.open(filename_sstr.str().c_str());

#ifdef AKANTU_USE_IOHELPER
  // without iohelper, this restart does not work and therefore the test is useless.
#else
  test_output << "Test does not make sense without iohelper" << std::endl;
  return EXIT_FAILURE;
#endif //AKANTU_USE_IOHELPER

  debug::setDebugLevel(dblWarning);
  akantu::initialize(argc, argv);

  /// load mesh
  Mesh my_mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("triangle_3.msh", my_mesh);

  /// build facet connectivity and surface id
  MeshUtils::buildFacets(my_mesh);
  MeshUtils::buildSurfaceID(my_mesh);

  UInt max_steps = 3;
  nb_nodes = my_mesh.getNbNodes();
  nb_elements = my_mesh.getNbElement(element_type);

  /// dump facet and surface information to paraview
#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperParaview dumper;
  dumper.SetMode(iohelper::TEXT);
  dumper.SetPoints(my_mesh.getNodes().storage(), dim, nb_nodes, "restart_triangle_3");
  dumper.SetConnectivity((int*)my_mesh.getConnectivity(element_type).storage(),
                         iohelper::TRIANGLE1, my_mesh.getNbElement(element_type), iohelper::C_MODE);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  /// declaration of model
  SolidMechanicsModel  my_model(my_mesh);
  /// model initialization
  my_model.initArrays();
  my_model.getForce().clear();
  my_model.getVelocity().clear();
  my_model.getAcceleration().clear();
  my_model.getDisplacement().clear();

  displacement = my_model.getDisplacement().storage();
  boundary = my_model.getBoundary().storage();

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

  master = 0;
  Surface impactor = 1;
  my_contact->addMasterSurface(master);
  my_contact->addImpactorSurfaceToMasterSurface(impactor, master);

  //FrictionCoefficient *fric_coef;
  //fric_coef = 
  new UniqueConstantFricCoef(*my_contact, master, 0.1);

  // contact information
  ContactRigid::SurfaceToImpactInfoMap::const_iterator it_imp_1;
  it_imp_1 = my_contact->getImpactorsInformation().find(master);
  ContactRigid::ImpactorInformationPerMaster * imp_info = it_imp_1->second;
  nb_components = imp_info->node_is_sticking->getNbComponent();

  my_model.updateCurrentPosition(); // neighbor structure uses current position for init
  my_contact->initNeighborStructure(master);
  my_contact->initSearch(); // does nothing so far

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(UInt s = 1; s <= max_steps; ++s) {

    std::cout << "passing step " << s << "/" << max_steps << "\r";
    std::cout.flush();

    /// apply a displacement to the slave body
    if(s == 2) {
      Real * coord = my_mesh.getNodes().storage();
      for(UInt n = 0; n < nb_nodes; ++n) {
        if(coord[n*dim + 0] > 0.5) {
          displacement[n*dim+0] = -0.02;
        }
      }
    }

    /// integration with contact and friction
    my_model.explicitPred();
    my_model.initializeUpdateResidualData();
    my_contact->solveContact();
    my_model.updateResidual(false);
    my_contact->avoidAdhesion();
    my_contact->frictionPredictor();
    my_model.updateAcceleration();
    my_model.explicitCorr();
    my_contact->frictionCorrector();
  }
  std::cout << std::endl;

  /// print contact info
  test_output << "CONTACT" << std::endl;
  printContact(my_contact);

  /// put restart information into a map
  std::map < std::string, ArrayBase* > output_map;
  my_contact->getRestartInformation(output_map, master);
  DumpRestart(my_model, output_map);

  test_output << "MAP" << std::endl;
  printRestartMap(output_map);
  //  freeMap(output_map);

  // ---------------- Reload contact from restart files ---------------------
  Contact * contact_restart = Contact::newContact(my_model,
                                                  _ct_rigid,
                                                  _cst_expli,
                                                  _cnst_regular_grid,
                                                  "contact_restart");

  ContactRigid * my_contact_restart = dynamic_cast<ContactRigid *>(contact_restart);
  my_contact_restart->initContact(false);
  my_contact_restart->addMasterSurface(master);
  my_contact_restart->addImpactorSurfaceToMasterSurface(impactor, master);
  std::map < std::string, ArrayBase* > input_map;
  loadRestartInformation(my_contact_restart, input_map);

  test_output << "RESTART MAP" << std::endl;
  printRestartMap(input_map);
  freeMap(input_map);

  test_output << "RESTART CONTACT" << std::endl;
  printContact(my_contact_restart);

  // finalize
  delete my_contact;
  delete my_contact_restart;
  finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
void restartReaderInit(iohelper::ReaderRestart & reader) {

#ifdef AKANTU_USE_IOHELPER
  reader.SetPoints("restart_test");
  reader.SetConnectivity(paraview_type);
  reader.AddNodeDataField("displacements");
  reader.AddNodeDataField("velocities");
  reader.AddNodeDataField("accelerations");
  reader.AddNodeDataField("forces");
  reader.AddNodeDataField("boundaries");
  reader.AddNodeDataField("active_impactor_nodes");
  reader.AddNodeDataField("master_element_type");
  reader.AddNodeDataField("master_normals");
  reader.AddNodeDataField("node_is_sticking");
  reader.AddNodeDataField("friction_forces");
  reader.AddNodeDataField("stick_positions");
  reader.AddNodeDataField("residual_forces");
  reader.AddNodeDataField("previous_velocities");
  reader.SetMode(iohelper::COMPRESSED);
  reader.SetPrefix("restart");
  reader.Init();
  reader.Read();

  // test if good number of node (-> good mesh)
  UInt reader_nb_nodes = reader.GetNumberNodes();
  if (reader_nb_nodes != nb_nodes)
    test_output << "WRONG MESH LOADED !!" << std::endl;

#endif //AKANTU_USE_IOHELPER
}

/* -------------------------------------------------------------------------- */
void loadRestartInformation(ContactRigid * contact, std::map < std::string, ArrayBase* > & map) {

#ifdef AKANTU_USE_IOHELPER
  // get the equilibrium state from the restart files
  iohelper::ReaderRestart * restart_reader = new iohelper::ReaderRestart();
  restartReaderInit(*restart_reader);
  memcpy(displacement,restart_reader->GetNodeDataField("displacements"),dim*nb_nodes*sizeof(Real));

  Real * tmp_r = restart_reader->GetNodeDataField("boundaries");
  Array<bool> * boundary_r = new Array<bool>(nb_nodes, dim, false);
  for (UInt i=0; i<nb_nodes; ++i) {
    for (UInt j=0; j<dim; ++j) {
      if(tmp_r[i*dim+j] > 0.01) {
        (*boundary_r)(i,j) = true;
      }
    }
  }
  memcpy(boundary, boundary_r->storage(), dim*nb_nodes*sizeof(bool));

  tmp_r = restart_reader->GetNodeDataField("active_impactor_nodes");
  Array<bool> * ai_nodes = new Array<bool>(nb_nodes, 1, false, "a_imp_nodes");
  for (UInt i=0; i<nb_nodes; ++i) {
    if(tmp_r[i] > 0.01)
      (*ai_nodes)(i) = true;
  }
  map["active_impactor_nodes"] = ai_nodes;

  tmp_r = restart_reader->GetNodeDataField("master_element_type");
  Array<ElementType> * et_nodes = new Array<ElementType>(nb_nodes, 1, _not_defined);
  for (UInt i=0; i<nb_nodes; ++i) {
    (*et_nodes)(i) = (ElementType)(tmp_r[i]);
  }
  map["master_element_type"] = et_nodes;

  tmp_r = restart_reader->GetNodeDataField("master_normals");
  Array<Real> * mn_nodes = new Array<Real>(0, dim);
  for (UInt i=0; i<nb_nodes; ++i) {
    Real normal[dim];
    for (UInt j=0; j<dim; ++j) {
      normal[j] = tmp_r[i*dim+j];
    }
    mn_nodes->push_back(normal);
  }
  map["master_normals"] = mn_nodes;

  tmp_r = restart_reader->GetNodeDataField("node_is_sticking");
  Array<bool> * is_nodes = new Array<bool>(nb_nodes, nb_components, false);
  for (UInt i=0; i<nb_nodes; ++i) {
    for (UInt j=0; j<nb_components; ++j) {
      if(tmp_r[i*nb_components+j] > 0.01)
        (*is_nodes)(i,j) = true;
    }
  }
  map["node_is_sticking"] = is_nodes;

  tmp_r = restart_reader->GetNodeDataField("friction_forces");
  Array<Real> * ff_nodes = new Array<Real>(0, dim);
  for (UInt i=0; i<nb_nodes; ++i) {
    Real normal[dim];
    for (UInt j=0; j<dim; ++j) {
      normal[j] = tmp_r[i*dim+j];
    }
    ff_nodes->push_back(normal);
  }
  map["friction_forces"] = ff_nodes;

  tmp_r = restart_reader->GetNodeDataField("stick_positions");
  Array<Real> * sp_nodes = new Array<Real>(0, dim);
  for (UInt i=0; i<nb_nodes; ++i) {
    Real position[dim];
    for (UInt j=0; j<dim; ++j) {
      position[j] = tmp_r[i*dim+j];
    }
    sp_nodes->push_back(position);
  }
  map["stick_positions"] = sp_nodes;

  tmp_r = restart_reader->GetNodeDataField("residual_forces");
  Array<Real> * rf_nodes = new Array<Real>(0, dim);
  for (UInt i=0; i<nb_nodes; ++i) {
    Real resid[dim];
    for (UInt j=0; j<dim; ++j) {
      resid[j] = tmp_r[i*dim+j];
    }
    rf_nodes->push_back(resid);
  }
  map["residual_forces"] = rf_nodes;

  tmp_r = restart_reader->GetNodeDataField("previous_velocities");
  Array<Real> * pv_nodes = new Array<Real>(0, dim);
  for (UInt i=0; i<nb_nodes; ++i) {
    Real pr_vel[dim];
    for (UInt j=0; j<dim; ++j) {
      pr_vel[j] = tmp_r[i*dim+j];
    }
    pv_nodes->push_back(pr_vel);
  }
  map["previous_velocities"] = pv_nodes;

  contact->setRestartInformation(map, master);
  delete boundary_r;
  // delete ai_nodes;
  // delete et_nodes;
  // delete mn_nodes;
  // delete is_nodes;
  // delete ff_nodes;
  // delete sp_nodes;
  // delete rf_nodes;
  // delete pv_nodes;
  delete restart_reader;
#endif //AKANTU_USE_IOHELPER
}

/* -------------------------------------------------------------------------- */
void DumpRestart(SolidMechanicsModel & my_model, std::map < std::string, ArrayBase* > & map) {
#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperRestart dumper;

  dumper.SetPoints(my_model.getFEM().getMesh().getNodes().storage(),
                   dim,nb_nodes,"restart_test");
  dumper.SetConnectivity((int *)my_model.getFEM().getMesh().getConnectivity(element_type).storage(),
                         paraview_type, nb_elements, iohelper::C_MODE);
  dumper.AddNodeDataField(my_model.getDisplacement().storage(),
                          dim, "displacements");
  dumper.AddNodeDataField(my_model.getVelocity().storage(),
                          dim, "velocities");
  dumper.AddNodeDataField(my_model.getAcceleration().storage()
                          ,dim, "accelerations");
  dumper.AddNodeDataField(my_model.getResidual().storage(),
                          dim, "forces");

  Array<Real> * boundary_r = new Array<Real>(nb_nodes, dim, 0.);
  for (UInt i=0; i<nb_nodes; ++i) {
    for (UInt j=0; j<dim; ++j) {
      if (boundary[i*dim+j]) {
        (*boundary_r)(i,j) = 1.;
      }
    }
  }
  dumper.AddNodeDataField(boundary_r->storage(), dim, "boundaries");

  // contact information
  std::map < std::string, ArrayBase* >::iterator it;

  it = map.find("active_impactor_nodes");
  Array<bool> * tmp_vec_b = (Array<bool> *)(it->second);
  Array<Real> * active_impactor_nodes = new Array<Real>(nb_nodes, 1, 0.);
  for (UInt i=0; i<nb_nodes; ++i) {
    if((*tmp_vec_b)(i))
      (*active_impactor_nodes)(i) = 1.;
  }
  dumper.AddNodeDataField(active_impactor_nodes->storage(), 1, "active_impactor_nodes");

  it = map.find("master_element_type");
  Array<ElementType> * tmp_vec_t = (Array<ElementType> *)(it->second);
  Array<Real> * master_element_type = new Array<Real>(nb_nodes, 1, (Real)_not_defined);
  for (UInt i=0; i<nb_nodes; ++i) {
    (*master_element_type)(i) = (Real)((*tmp_vec_t)(i));
  }
  dumper.AddNodeDataField(master_element_type->storage(), 1, "master_element_type");

  it = map.find("master_normals");
  dumper.AddNodeDataField(((Array<Real> *)it->second)->storage(), dim, "master_normals");

  it = map.find("node_is_sticking");
  tmp_vec_b = (Array<bool> *)(it->second);
  Array<Real> * node_is_sticking = new Array<Real>(nb_nodes, 2, 0.);
  for (UInt i=0; i<nb_nodes; ++i) {
    for (UInt j=0; j<2; ++j) {
      if((*tmp_vec_b)(i,j))
        (*node_is_sticking)(i,j) = 1.;
    }
  }
  dumper.AddNodeDataField(node_is_sticking->storage(), 2, "node_is_sticking");

  it = map.find("friction_forces");
  dumper.AddNodeDataField(((Array<Real> *)it->second)->storage(), dim, "friction_forces");

  it = map.find("stick_positions");
  dumper.AddNodeDataField(((Array<Real> *)it->second)->storage(), dim, "stick_positions");

  it = map.find("residual_forces");
  dumper.AddNodeDataField(((Array<Real> *)it->second)->storage(), dim, "residual_forces");
  it = map.find("previous_velocities");
  dumper.AddNodeDataField(((Array<Real> *)it->second)->storage(), dim, "previous_velocities");

  dumper.SetMode(iohelper::COMPRESSED);
  dumper.SetPrefix("restart");
  dumper.Init();

  dumper.Dump();

  delete boundary_r;
  delete active_impactor_nodes;
  delete master_element_type;
  delete node_is_sticking;
#endif //AKANTU_USE_IOHELPER
}

/* -------------------------------------------------------------------------- */
void printRestartMap(std::map < std::string, ArrayBase* > & map) {
  /// access all vectors in the map
  std::map < std::string, ArrayBase* >::iterator it;

  it = map.find("active_impactor_nodes");
  Array<bool> * ai_nodes = dynamic_cast<Array<bool> *>(it->second);
  if(it == map.end()) {
    test_output << "could not find map entry for active impactor nodes" << std::endl;
  }

  it = map.find("master_element_type");
  Array<ElementType> * et_nodes = dynamic_cast<Array<ElementType> *>(it->second);
  if(it == map.end()) {
    test_output << "could not find map entry master element type" << std::endl;
  }

  it = map.find("master_normals");
  Array<Real> * mn_nodes = dynamic_cast<Array<Real> *>(it->second);
  if(it == map.end()) {
    test_output << "could not find map entry for master normals" << std::endl;
  }

  it = map.find("node_is_sticking");
  Array<bool> * is_nodes = dynamic_cast<Array<bool> *>(it->second);
  if(it == map.end()) {
    test_output << "could not find map entry node is sticking" << std::endl;
  }

  it = map.find("friction_forces");
  Array<Real> * ff_nodes = dynamic_cast<Array<Real> *>(it->second);
  if(it == map.end()) {
    test_output << "could not find map entry friction forces" << std::endl;
  }

  it = map.find("stick_positions");
  Array<Real> * sp_nodes = dynamic_cast<Array<Real> *>(it->second);
  if(it == map.end()) {
    test_output << "could not find map entry stick positions" << std::endl;
  }

  it = map.find("residual_forces");
  Array<Real> * rf_nodes = dynamic_cast<Array<Real> *>(it->second);
  if(it == map.end()) {
    test_output << "could not find map entry residual forces" << std::endl;
  }

  it = map.find("previous_velocities");
  Array<Real> * pv_nodes = dynamic_cast<Array<Real> *>(it->second);
  if(it == map.end()) {
    test_output << "could not find map entry previous velocities" << std::endl;
  }

  /// print all information in restart map
  test_output << "Active impactor nodes (map):" << std::endl << std::endl;
  for (UInt i=0; i<nb_nodes; ++i) {
    if ((*ai_nodes)(i)) {
      test_output << "node: " << i << ", master element type: " << (*et_nodes)(i) << std::endl;
      for (UInt d=0; d<dim; ++d) {
        test_output << "Direction " << d << std::endl;
        test_output << "  master normal = " << (*mn_nodes)(i,d) << std::endl;
        test_output << "  friction force = " << (*ff_nodes)(i,d) << std::endl;
        test_output << "  stick position = " << (*sp_nodes)(i,d) << std::endl;
        test_output << "  residual force = " << (*rf_nodes)(i,d) << std::endl;
        test_output << "  previous velocity = " << (*pv_nodes)(i,d) << std::endl;
      }
      for (UInt j=0; j<2; ++j) {
        test_output << "stick " << j << ": " << (*is_nodes)(i,j) << std::endl;
      }
      test_output << std::endl;
    }
  }
}

/* -------------------------------------------------------------------------- */
void printContact(ContactRigid * contact) {
  ContactRigid::SurfaceToImpactInfoMap::const_iterator it_imp;
  it_imp = contact->getImpactorsInformation().find(master);
  ContactRigid::ImpactorInformationPerMaster * impactor_info = it_imp->second;

  UInt * active_nodes = impactor_info->active_impactor_nodes->storage();
  ElementType * element_type_imp = &(*impactor_info->master_element_type)[0];
  Real * master_normal = impactor_info->master_normals->storage();
  bool * node_stick = impactor_info->node_is_sticking->storage();
  Real * friction_force = impactor_info->friction_forces->storage();
  Real * stick_position = impactor_info->stick_positions->storage();
  Real * residual_force = impactor_info->residual_forces->storage();
  Real * previous_velocity = impactor_info->previous_velocities->storage();

  UInt nb_active_nodes = impactor_info->active_impactor_nodes->getSize();
  test_output << "Active impactor nodes (contact):" << std::endl << std::endl;
  for (UInt i=0; i<nb_active_nodes; ++i, ++active_nodes, ++element_type_imp) {
    test_output << "node: " <<  (*active_nodes)
                << ", master element type: " << *element_type_imp << std::endl;
    for (UInt d=0; d<dim; ++d, ++master_normal,
           ++friction_force,
           ++stick_position,
           ++residual_force,
           ++previous_velocity) {
      test_output << "Direction " << d << std::endl;
      test_output << "  master normal = " << *master_normal << std::endl;
      test_output << "  friction force = " << *friction_force << std::endl;
      test_output << "  stick position = " << *stick_position << std::endl;
      test_output << "  residual force = " << *residual_force << std::endl;
      test_output << "  previous velocity = " << *previous_velocity << std::endl;
    }
    for (UInt j=0; j<2; ++j, ++node_stick) {
      test_output << "stick " << j << ": " << *node_stick << std::endl;
    }
    test_output << std::endl;
  }
}


/* -------------------------------------------------------------------------- */
void freeMap(std::map < std::string, ArrayBase* > & map) {
  std::map < std::string, ArrayBase* >::iterator it = map.begin();
  std::map < std::string, ArrayBase* >::iterator end = map.end();

  for (; it != end; ++it) {
    delete it->second;
  }
}
