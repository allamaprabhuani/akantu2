/**
 * @file   test_contact_2d_neighbor_structure.cc
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @date   Thu Dec  9 10:07:58 2010
 *
 * @brief  Test neighbor structure for 2d with linear triangles
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
#include "contact_2d_explicit.hh"
#include "contact_neighbor_structure.hh"

using namespace akantu;

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"

static void initParaview(Mesh & mesh);
static void initParaviewSurface(Mesh & mesh);
static void printParaviewSurface(Mesh & mesh, const NeighborList & my_neighbor_list);
double * facet_id;
double * node_id;
iohelper::DumperParaview dumper_surface;
#endif //AKANTU_USE_IOHELPER



int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);

  int spatial_dimension = 2;
  Real time_factor = 0.2;

  /// load mesh
  Mesh mesh(spatial_dimension);
  MeshIOMSH mesh_io;
  mesh_io.read("squares.msh", mesh);

  SolidMechanicsModel * model = new SolidMechanicsModel(mesh);

  UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  UInt nb_elements = model->getFEM().getMesh().getNbElement(_triangle_3);

  /// model initialization
  model->initModel();
  model->initVectors();

  model->readMaterials("materials.dat");
  model->initMaterials();

  model->initExplicit();
  std::cout << model->getMaterial(0) << std::endl;

  model->assembleMassLumped();

  /// set vectors to zero
  memset(model->getForce().values,        0,
	 spatial_dimension*nb_nodes*sizeof(Real));
  memset(model->getVelocity().values,     0,
	 spatial_dimension*nb_nodes*sizeof(Real));
  memset(model->getAcceleration().values, 0,
	 spatial_dimension*nb_nodes*sizeof(Real));
  memset(model->getDisplacement().values, 0,
	 spatial_dimension*nb_nodes*sizeof(Real));
  memset(model->getResidual().values, 0,
	 spatial_dimension*nb_nodes*sizeof(Real));
  memset(model->getMaterial(0).getStrain(_triangle_3).values, 0,
	 spatial_dimension*spatial_dimension*nb_elements*sizeof(Real));
  memset(model->getMaterial(0).getStress(_triangle_3).values, 0,
	 spatial_dimension*spatial_dimension*nb_elements*sizeof(Real));

  /// Paraview Helper
#ifdef AKANTU_USE_IOHELPER
  // initParaview(*model);
#endif //AKANTU_USE_IOHELPER

//   /// dump surface information to paraview
// #ifdef AKANTU_USE_IOHELPER
//   iohelper::DumperParaview dumper;
//   dumper.SetMode(iohelper::TEXT);
  
//   dumper.SetPoints(mesh.getNodes().values, spatial_dimension, nb_nodes, "triangle_3_test-surface-extraction");
//   dumper.SetConnectivity((int*)mesh.getConnectivity(_triangle_3).values,
//    			 iohelper::TRIANGLE1, mesh.getNbElement(_triangle_3), iohelper::C_MODE);
//   dumper.SetPrefix("paraview/");
//   dumper.Init();
//   dumper.Dump();
// #endif //AKANTU_USE_IOHELPER

  Real time_step = model->getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);
  

  /// contact declaration
  Contact * contact = Contact::newContact(*model, 
					  _ct_2d_expli, 
					  _cst_2d_expli, 
					  _cnst_2d_grid);

  Contact2dExplicit * my_contact = dynamic_cast<Contact2dExplicit *>(contact);

  my_contact->initContact(true);
  my_contact->setFrictionCoefficient(0.);
  my_contact->initNeighborStructure();

  /// get master surfaces with associated neighbor list
  //  const std::vector<Surface> & master_surfaces = my_contact->getMasterSurfaces();
  std::vector<Surface>::iterator it;
  // for (it = master_surfaces.begin(); it != master_surfaces.end(); ++it) {

#ifdef AKANTU_USE_IOHELPER
  initParaview(mesh);
  initParaviewSurface(mesh);
#endif //AKANTU_USE_IOHELPER
  

  UInt nb_surfaces = mesh.getNbSurfaces();
  for (UInt s = 0; s < nb_surfaces; ++s) {

    const NeighborList & my_neighbor_list = my_contact->getContactSearch().getContactNeighborStructure(s).getNeighborList();

#ifdef AKANTU_USE_IOHELPER
    printParaviewSurface(mesh, my_neighbor_list);
#endif //AKANTU_USE_IOHELPER
  
    UInt nb_impactors = my_neighbor_list.impactor_nodes.getSize();
    UInt * impactors_val = my_neighbor_list.impactor_nodes.values;
  
    UInt * node_facet_off_val = my_neighbor_list.facets_offset(_segment_2).storage();
    UInt * node_facet_val = my_neighbor_list.facets(_segment_2).storage();

    /// print impactor nodes
    std::cout << "Master surface " << s << " has " <<  nb_impactors << " impactor nodes:" << std::endl;
    for (UInt i = 0; i < nb_impactors; ++i) {
      std::cout << " node " << impactors_val[i] << " : ";
      for (UInt j = node_facet_off_val[i]; j < node_facet_off_val[i+1]; ++j)
	std::cout << node_facet_val[j] << " ";
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  // }


#ifdef AKANTU_USE_IOHELPER
  delete [] node_id;
  delete [] facet_id;
#endif //AKANTU_USE_IOHELPER
  
  delete my_contact;
  delete model;
  finalize();

  return EXIT_SUCCESS;
}


/// Paraview prints
#ifdef AKANTU_USE_IOHELPER

static void initParaview(Mesh & mesh) {

  iohelper::DumperParaview dumper;
  dumper.SetMode(iohelper::TEXT);

  UInt  nb_nodes = mesh.getNbNodes();
  dumper.SetPoints(mesh.getNodes().values, 2, nb_nodes, "test-2d-neighbor");
  dumper.SetConnectivity((int*)mesh.getConnectivity(_triangle_3).values,
   			 iohelper::TRIANGLE1, mesh.getNbElement(_triangle_3), iohelper::C_MODE);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
}


static void initParaviewSurface(Mesh & mesh) {
  
  iohelper::DumperParaview dumper_surface;

  dumper_surface.SetMode(iohelper::TEXT);

  UInt  nb_nodes = mesh.getNbNodes();
  dumper_surface.SetPoints(mesh.getNodes().values, 2, nb_nodes, "test-2d-neighbor-surface");


  dumper_surface.SetConnectivity((int *)mesh.getConnectivity(_segment_2).values,
			       iohelper::LINE1, mesh.getNbElement(_segment_2), iohelper::C_MODE);

  
  facet_id = new double [mesh.getConnectivity(_segment_2).getSize()];
  memset(facet_id, 0, mesh.getConnectivity(_segment_2).getSize()*sizeof(double));

  node_id = new double [nb_nodes];
  memset(node_id, 0, nb_nodes*sizeof(double));

  dumper_surface.AddElemDataField(facet_id, 1, "master_segments");
  dumper_surface.AddNodeDataField(node_id, 1, "slave_nodes");

  dumper_surface.SetPrefix("paraview/");
  dumper_surface.Init();
  dumper_surface.Dump();
  // delete [] facet_id;
  // delete [] node_id;
  // return  dumper_surface;
}

static void printParaviewSurface(Mesh & mesh, const NeighborList & my_neighbor_list) {

  UInt nb_impactors = my_neighbor_list.impactor_nodes.getSize();
  UInt * impactors_val = my_neighbor_list.impactor_nodes.values;
  
  UInt nb_facets = my_neighbor_list.facets(_segment_2).getSize();
  UInt * node_facet_val = my_neighbor_list.facets(_segment_2).storage();

  // double * node_id = new double [mesh.getNbNodes()];
  memset(node_id, 0, mesh.getNbNodes()*sizeof(double));

  // double * facet_id = new double [mesh.getConnectivity(_segment_2).getSize()];
  memset(facet_id, 0, mesh.getConnectivity(_segment_2).getSize()*sizeof(double));

  for (UInt n = 0; n < nb_impactors; ++n)
    node_id[impactors_val[n]] = 1.;

  for (UInt el = 0; el < nb_facets; ++el)
    facet_id[node_facet_val[el]] = 1.;

  dumper_surface.Dump();
}

#endif //AKANTU_USE_IOHELPER
