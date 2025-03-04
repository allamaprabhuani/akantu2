/**
 * Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
#include "material_cohesive.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>
#include <limits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {

  initialize("input_file.dat", argc, argv);

  Math::setTolerance(1e-15);

  const Int spatial_dimension = 3;

  Mesh mesh(spatial_dimension);

  mesh.read("3d_spherical_inclusion.msh");

  SolidMechanicsModelCohesive model(mesh);

  auto && material_selector =
      std::make_shared<MeshDataMaterialCohesiveSelector>(model);
  material_selector->setFallback(model.getMaterialSelector());
  model.setMaterialSelector(material_selector);

  model.initFull();

  std::vector<std::string> surfaces_name = {"interface", "coh1", "coh2",
                                            "coh3",      "coh4", "coh5"};
  Int nb_surf = surfaces_name.size();

  for (const auto & type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_cohesive)) {
    for (Int i = 0; i < nb_surf; ++i) {

      auto expected_insertion = mesh.getElementGroup(surfaces_name[i])
                                    .getElements(mesh.getFacetType(type))
                                    .size();
      const auto & mat = model.getMaterial(surfaces_name[i]);
      auto inserted_elements = mat.getElementFilter()(type).size();
      if (not(expected_insertion == inserted_elements)) {
        std::cout << "!!! Mismatch in insertion of surface named "
                  << surfaces_name[i] << " --> " << inserted_elements
                  << " inserted elements out of " << expected_insertion
                  << std::endl;
        return 1;
      }
    }
  }

  return 0;
}
