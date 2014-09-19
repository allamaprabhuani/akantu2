/**
 * @file   acurnier_3D_2.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date   Thu Sep 10 16:30:00 2014
 *
 * @brief  This file implements the fifth simple test suggested by Alain
 * Curnier for the verification of the implicit contact implementation of the 
 * Augmented lagrangian formulation
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

#include "contact_impl.hh"

using namespace akantu;


int main(int argc, char *argv[]) {
  
  // set dimension
  static const UInt dim = 3;
  
  // type definitions
  typedef SolidMechanicsModel model_type;
  typedef Contact<dim, MasterAssignator,
  SelectResolution<_static, _augmented_lagrangian> >
  contact_type;
  
  initialize("material.dat", argc, argv);
  
  // create meshes
  Mesh mesh(dim);
  
  // read meshes
  mesh.read("acurnier_3D_2.msh");
  
  // create models
  model_type model(mesh);
  
  SolidMechanicsModelOptions opt(_static);
  
  MeshDataMaterialSelector<std::string> material_selector("physical_names", model);
  model.setMaterialSelector(material_selector);


  // initialize material
  model.initFull(opt);
  model.updateCurrentPosition();

  // create data structure that holds contact data
  contact_type cd(argc, argv, model);
  
  // set Paraview output resluts
  model.setBaseName("contact");
  model.addDumpFieldVector("displacement");
  
  mesh.createGroupsFromMeshData<std::string>("physical_names");
  
  cd.addSlave(5);
  cd.addSlave(6);
  cd.addSlave(11);
  cd.addArea(5, 1.0);
  cd.addArea(6, 1.0);
  cd.addArea(11, 1.0);
  
  // add master surface to find pairs
  cd.searchSurface("Contact");
  
  model.applyBC(BC::Dirichlet::FixedValue(0., _x), "Top");
  model.applyBC(BC::Dirichlet::FixedValue(0., _z), "Top");

  // fix entire contact body
  model.applyBC(BC::Dirichlet::FixedValue(0., _x), "Bottom");
  model.applyBC(BC::Dirichlet::FixedValue(0., _y), "Bottom");
  model.applyBC(BC::Dirichlet::FixedValue(0., _z), "Bottom");
  
  Real U = 0.5;
  Real Du = 0.01;
  for (Real u = Du; u <= U; u += Du) {
    
    model.applyBC(BC::Dirichlet::FixedValue(-u, _y), "Top");
    
    // solve contact step (no need to call solve on the model object)
    solveContactStep<_generalized_newton>(cd);
  }
  
  cout<<"Force: "<< cd.getForce()<<endl;
  
  // finalize simulation
  finalize();
  return EXIT_SUCCESS;
}

