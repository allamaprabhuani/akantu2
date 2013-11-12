/**
 * @file   test_cohesive_intrinsic_tetrahedron.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Tue Aug 20 14:37:18 2013
 *
 * @brief  Test for cohesive elements
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
#include "mesh_utils.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "material.hh"
#include "dumper_paraview.hh"

/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  initialize(argc, argv);

  //  debug::setDebugLevel(dblDump);
  ElementType type = _tetrahedron_10;

  const UInt spatial_dimension = 3;
  const UInt max_steps = 100;

  Mesh mesh(spatial_dimension);
  mesh.read("tetrahedron_full.msh");


  /* ------------------------------------------------------------------------ */
  /* Facet part                                                               */
  /* ------------------------------------------------------------------------ */

  // Array<Real> limits(spatial_dimension, 2);
  // limits(0, 0) = -100;
  // limits(0, 1) = 100;
  // limits(1, 0) = -100;
  // limits(1, 1) = 100;
  // limits(2, 0) = -100;
  // limits(2, 1) = 100;

  // MeshUtils::insertIntrinsicCohesiveElementsInArea(mesh, limits);

  Mesh mesh_facets(spatial_dimension, mesh.getNodes(), "mesh_facets");
  MeshUtils::buildAllFacets(mesh, mesh_facets);

  MeshUtils::resetFacetToDouble(mesh_facets);

  ByElementTypeArray<bool> facet_insertion("facet_insertion", "");
  mesh_facets.initByElementTypeArray(facet_insertion, 1, spatial_dimension - 1,
				     false, _ek_regular, true);

  const ElementType type_facet = Mesh::getFacetType(type);
  Array<bool> & f_insertion = facet_insertion(type_facet);
  Array<std::vector<Element> > & element_to_facet
    = mesh_facets.getElementToSubelement(type_facet);

  UInt nb_facet = mesh_facets.getNbElement(type_facet);

  for (UInt f = 0; f < nb_facet; ++f) {

    if (element_to_facet(f)[1] == ElementNull) continue;

    f_insertion(f) = true;

    MeshUtils::insertCohesiveElements(mesh,
				      mesh_facets,
				      facet_insertion,
				      false);
  }

  /* ------------------------------------------------------------------------ */
  /* End of facet part                                                        */
  /* ------------------------------------------------------------------------ */

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull("material.dat");
  Real time_step = model.getStableTimeStep()*0.8;
  model.setTimeStep(time_step);
  //  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();

  model.updateResidual();

  model.setBaseName("intrinsic_tetrahedron_fragmentation");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("residual"    );
  model.addDumpField("stress");
  model.addDumpField("strain");
  model.dump();

  DumperParaview dumper("cohesive_elements_tetrahedron_fragmentation");
  dumper.registerMesh(mesh, spatial_dimension, _not_ghost, _ek_cohesive);
  DumperIOHelper::Field * cohesive_displacement =
    new DumperIOHelper::NodalField<Real>(model.getDisplacement());
  cohesive_displacement->setPadding(3);
  dumper.registerField("displacement", cohesive_displacement);
  // dumper.registerField("damage", new DumperIOHelper::
  // 		       HomogenizedField<Real,
  // 					DumperIOHelper::InternalMaterialField>(model,
  // 									       "damage",
  // 									       spatial_dimension,
  // 									       _not_ghost,
  // 									       _ek_cohesive));

  dumper.dump();

  /// update displacement
  UInt nb_element = mesh.getNbElement(type);
  UInt nb_nodes = mesh.getNbNodes();
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  Real * bary = new Real[spatial_dimension];

  const Array<UInt> & connectivity = mesh.getConnectivity(type);
  Array<Real> & displacement = model.getDisplacement();
  Array<bool> update(nb_nodes);


  for (UInt s = 0; s < max_steps; ++s) {
    Real increment = s / 10.;
    update.clear();

    for (UInt el = 0; el < nb_element; ++el) {
      mesh.getBarycenter(el, type, bary);
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	UInt node = connectivity(el, n);
	if (!update(node)) {
	  for (UInt dim = 0; dim < spatial_dimension; ++dim) {
	    displacement(node, dim) = increment * bary[dim];
	    update(node) = true;
	  }
	}
      }
    }

    if (s % 10 == 0) {
      model.dump();
      dumper.dump();
    }
  }

  delete[] bary;

  if (nb_nodes != nb_element * Mesh::getNbNodesPerElement(type)) {
    std::cout << "Wrong number of nodes" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  finalize();
  std::cout << "OK: test_cohesive_intrinsic_tetrahedron was passed!" << std::endl;
  return EXIT_SUCCESS;
}
