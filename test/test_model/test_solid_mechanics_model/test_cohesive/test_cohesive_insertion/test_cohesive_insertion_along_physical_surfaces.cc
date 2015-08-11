/**
 * @file   test_cohesive_insertion_along_physical_surfaces.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Fri Aug  7 09:07:44 2015
 *
 * @brief Test intrinsic insertion of cohesive elements along physical surfaces 
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
#include <limits>
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "material.hh"
#include "material_cohesive.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[]) {

  initialize("input_file.dat", argc, argv);

  Math::setTolerance(1e-15);
  
  const UInt spatial_dimension = 3;

  Mesh mesh(spatial_dimension);

  mesh.read("3d_spherical_inclusion.msh");
  
  SolidMechanicsModelCohesive model(mesh);
  mesh.createGroupsFromMeshData<std::string>("physical_names");
  model.initFull(SolidMechanicsModelCohesiveOptions(_static));
  
  std::vector<std::string> surfaces_name = {"interface", "coh1", "coh2", "coh3", "coh4", "coh5"};
  UInt nb_surf = surfaces_name.size();

  Mesh::type_iterator it  = mesh.firstType(spatial_dimension, _not_ghost, _ek_cohesive);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, _not_ghost, _ek_cohesive);

  for(; it != end; ++it) {
    for (UInt i = 0; i < nb_surf; ++i) {
      
      UInt expected_insertion = mesh.getElementGroup(surfaces_name[i]).getElements(mesh.getFacetType(*it)).getSize();
      UInt inserted_elements = model.getMaterial(surfaces_name[i]).getElementFilter()(*it).getSize();
      AKANTU_DEBUG_ASSERT((expected_insertion == inserted_elements),
			  std::endl << "!!! Mismatch in insertion of surface named " 
			  << surfaces_name[i] 
			  << " --> " << inserted_elements << " inserted elements out of " 
			  << expected_insertion << std::endl);
    }      
  }

  /*std::string paraview_folder = "paraview/test_intrinsic_insertion_along_physical_surfaces/"; 
  model.setDirectory(paraview_folder);
  model.setBaseName("bulk");
  model.addDumpField("partitions");
  model.dump();
  model.setDirectoryToDumper("cohesive elements", paraview_folder);
  model.setBaseNameToDumper("cohesive elements", "one_cohesive_element");
  model.addDumpFieldToDumper("cohesive elements", "partitions");
  model.addDumpFieldToDumper("cohesive elements", "material_index");
  model.dump("cohesive elements");
  */

  model.assembleStiffnessMatrix();

  finalize();

  return EXIT_SUCCESS;
}
