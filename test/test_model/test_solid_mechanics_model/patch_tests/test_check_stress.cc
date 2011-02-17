/**
 * @file   test_check_stress.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Wed Feb 16 13:56:42 2011
 *
 * @brief  patch test for elastic material in solid mechanics model
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

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  int dim = 2;
  const ElementType element_type = _triangle_3;
  const UInt paraview_type = TRIANGLE1;
  
  UInt imposing_steps = 1000;
  Real max_displacement = -0.01;

  UInt damping_steps = 300000;
  UInt damping_interval = 50;
  Real damping_ratio = 0.99;

  UInt additional_steps = 20000;
  UInt max_steps = damping_steps + additional_steps;
  std::cout << "The number of time steps is: " << max_steps << std::endl; 

  /// load mesh
  Mesh my_mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("square_check_stress.msh", my_mesh);

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
  my_model.readMaterials("material_check_stress.dat");
  my_model.initMaterials();

  UInt nb_element = my_model.getFEM().getMesh().getNbElement(element_type);

  Real time_step = my_model.getStableTimeStep();
  my_model.setTimeStep(time_step/10.);

  my_model.assembleMassLumped();

  // boundary conditions
  Vector<UInt> * top_nodes  = new Vector<UInt>(0, 1);
  Vector<UInt> * base_nodes = new Vector<UInt>(0, 1);
  Real * coordinates = my_mesh.getNodes().values;
  Real * displacement = my_model.getDisplacement().values;
  bool * boundary = my_model.getBoundary().values;
  for(UInt n = 0; n < nb_nodes; ++n) {
    Real y_coord = coordinates[n*dim + 1];
    if (y_coord > 0.99) {
      boundary[n*dim + 1] = true;
      top_nodes->push_back(n);
    }
    else if (y_coord < 0.01) {
      boundary[n*dim + 1] = true;
      base_nodes->push_back(n);
    }
  }
  UInt * top_nodes_val  = top_nodes->values;
  UInt * base_nodes_val = base_nodes->values;
  
  Real * velocity_val = my_model.getVelocity().values;

  /*
#ifdef AKANTU_USE_IOHELPER
  /// initialize the paraview output
  DumperParaview dumper;
  dumper.SetMode(TEXT);
  dumper.SetPoints(my_model.getFEM().getMesh().getNodes().values, dim, nb_nodes, "coord_square_check_stress");
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
  dumper.SetPrefix("paraview/check_stress_test/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  */

  /*
  std::ofstream out_info;
  out_info.open("check_force.csv");
  out_info << "%id,ftop,fcont" << std::endl;

  std::ofstream energy;
  energy.open("check_force_energy.csv");
  energy << "%id,kin,pot,tot" << std::endl;
  */

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(UInt s = 1; s <= max_steps; ++s) {

    if(s % 10000 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;

    // impose normal displacement
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
   
    my_model.updateResidual();

    /*
    // find the total force applied at the imposed displacement surface (top)
    Real * residual = my_model.getResidual().values; 
    Real top_force = 0.;
    Real base_force = 0.;
    for(UInt n=0; n<top_nodes->getSize(); ++n) {
      UInt node = top_nodes_val[n];
      top_force += residual[node*dim + 1];
    }
    for(UInt n=0; n<base_nodes->getSize(); ++n) {
      UInt node = base_nodes_val[n];
      base_force += residual[node*dim + 1];
    }
    out_info << s << "," << top_force << "," << base_force << std::endl;
    */

    my_model.updateAcceleration();
    my_model.explicitCorr();

    /*
    Real epot = my_model.getPotentialEnergy();
    Real ekin = my_model.getKineticEnergy();
    energy << s << "," << ekin << "," << epot << "," << ekin+epot << std::endl;
    */
    /*
#ifdef AKANTU_USE_IOHELPER
    if(s % 10000 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
    */
  }

  /*
  out_info.close();
  energy.close();
  */
  
  UInt check_element = 0;
  UInt quadrature_point = 0;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element_type);
  UInt nb_quadrature_points = my_model.getFEM().getNbQuadraturePoints(element_type);
  Real * strains  = my_model.getMaterial(0).getStrain(element_type).values;
  Real * stresses = my_model.getMaterial(0).getStress(element_type).values;
  Real * strain_val = &strains[check_element*quadrature_point*dim*dim + quadrature_point];
  Real * stress_val = &stresses[check_element*quadrature_point*dim*dim + quadrature_point];
  
  if(abs(strain_val[0] - 1.004285714) > 1e-9) {
    std::cout << "strain[0] = " << strain_val[0] << " but should be = 1.004285714" << std::endl;
    return EXIT_FAILURE;
  }

  if(abs(strain_val[1]) > 1e-15) {
    std::cout << "strain[1] = " << strain_val[1] << " but should be = 0" << std::endl;
    return EXIT_FAILURE;
  }

  if(abs(strain_val[2]) > 1e-15) {
    std::cout << "strain[2] = " << strain_val[2] << " but should be = 0" << std::endl;
    return EXIT_FAILURE;
  }

  if(abs(strain_val[3] - 0.99) > 1e-9) {
    std::cout << "strain[3] = " << strain_val[3] << " but should be = 0.99" << std::endl;
    return EXIT_FAILURE;
  }

  if(abs(stress_val[0]) > 1e-4) {
    std::cout << "stress[0] = " << stress_val[0] << " but should be = 0" << std::endl;
    return EXIT_FAILURE;
  }

  if(abs(stress_val[1]) > 1e-4) {
    std::cout << "stress[1] = " << stress_val[1] << " but should be = 0" << std::endl;
    return EXIT_FAILURE;
  }

  if(abs(stress_val[2]) > 1e-4) {
    std::cout << "stress[2] = " << stress_val[2] << " but should be = 0" << std::endl;
    return EXIT_FAILURE;
  }

  if(abs(stress_val[3] + 2.30769e9) > 1e5) {
    std::cout << "stress[3] = " << stress_val[3] << " but should be = -2.30769e9" << std::endl;
    return EXIT_FAILURE;
  }

  finalize();

  return EXIT_SUCCESS;
}
