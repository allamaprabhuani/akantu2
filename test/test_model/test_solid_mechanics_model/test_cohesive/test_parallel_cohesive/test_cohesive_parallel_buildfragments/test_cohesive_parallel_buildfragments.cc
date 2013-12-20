/**
 * @file   test_cohesive_parallel_buildfragments.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Tue Dec 17 12:03:03 2013
 *
 * @brief  Test to build fragments in parallel
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

/* -------------------------------------------------------------------------- */

#include <limits>
#include <fstream>
#include <iostream>


/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "material_cohesive.hh"

//#include "io_helper.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

void verticalInsertionLimit(SolidMechanicsModelCohesive &);
void displaceElements(SolidMechanicsModelCohesive &, const Real, const Real);

int main(int argc, char *argv[]) {
  initialize(argc, argv);

  debug::setDebugLevel(dblWarning);

  const UInt spatial_dimension = 2;
  const UInt total_nb_fragment = 2;

  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    // Read the mesh
    mesh.read("mesh.msh");

    /// partition the mesh
    MeshUtils::purifyMesh(mesh);
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }

  SolidMechanicsModelCohesive model(mesh);
  model.initParallel(partition, NULL, true);

  delete partition;

  /// model initialization
  model.initFull("material.dat",
		 SolidMechanicsModelCohesiveOptions(_explicit_lumped_mass, true));

  mesh.computeBoundingBox();
  Real L = mesh.getXMax() - mesh.getXMin();
  Real h = mesh.getYMax() - mesh.getYMin();
  Real rho = model.getMaterial("bulk").getParam<Real>("rho");

  Real theoretical_mass = L * h * rho;
  Real frag_theo_mass = theoretical_mass / total_nb_fragment;

  model.setBaseName("extrinsic");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("residual"    );
  model.addDumpField("stress");
  model.addDumpField("strain");
  model.addDumpField("partitions");
  model.dump();

  /// set check facets
  verticalInsertionLimit(model);

  model.assembleMassLumped();
  model.synchronizeBoundaries();

  /// impose initial displacement
  Array<Real> & displacement = model.getDisplacement();
  const Array<Real> & position = mesh.getNodes();
  UInt nb_nodes = mesh.getNbNodes();

  for (UInt n = 0; n < nb_nodes; ++n)
    displacement(n, 0) = position(n, 0) * 10;

  model.updateResidual();
  model.checkCohesiveStress();
  model.dump();

  {
    ElementType cohesive_type = _cohesive_2d_4;
    UInt nb_element = mesh.getNbElement(cohesive_type);
    UInt nb_element_ghost = mesh.getNbElement(cohesive_type, _ghost);

    std::cout << "Proc " << prank << ": " << nb_element
  	      << " _not_ghost, " << nb_element_ghost
  	      << " _ghost" << std::endl;

    // Vector<Real> barycenter(spatial_dimension);

    // for (UInt el = 0; el < nb_element; ++el) {
    //   mesh.getBarycenter(el, cohesive_type, barycenter.storage());

    //   std::cout << "  " << el << ": (" << barycenter(0)
    // 		<< ", " << barycenter(1) << ")" << std::endl;
    // }

    // std::cout << "Ghost:" << std::endl;

    // for (UInt el = 0; el < nb_element_ghost; ++el) {
    //   mesh.getBarycenter(el, cohesive_type, barycenter.storage(), _ghost);

    //   std::cout << "  " << el << ": (" << barycenter(0)
    // 		<< ", " << barycenter(1) << ")" << std::endl;
    // }
  }

  const Array<Real> & fragment_mass = model.getFragmentsMass();

  Real el_size = L / total_nb_fragment;
  Real lim = -L/2 + el_size * 0.75;

  /// displace one fragment each time
  for (UInt frag = 1; frag <= total_nb_fragment; ++frag) {
    if (prank == 0)
      std::cout << "Generating fragment: " << frag << std::endl;

    model.computeFragmentsData();

    /// check number of fragments
    UInt nb_fragment_num = model.getNbFragment();

    if (nb_fragment_num != frag) {
      AKANTU_DEBUG_ERROR("The number of fragments is wrong! Numerical: " << nb_fragment_num << " Theoretical: " << frag);
      return EXIT_FAILURE;
    }

    /// check mass computation
    if (frag < total_nb_fragment) {
      Real total_mass = 0.;
      UInt small_fragments = 0;

      for (UInt f = 0; f < nb_fragment_num; ++f) {
	if (Math::are_float_equal(fragment_mass(f, 0), frag_theo_mass)) {
	  ++small_fragments;
	  total_mass += frag_theo_mass;
	}
      }

      if (small_fragments != nb_fragment_num - 1) {
	AKANTU_DEBUG_ERROR("The number of small fragments is wrong!");
	return EXIT_FAILURE;
      }

      if (!Math::are_float_equal(total_mass,
				 small_fragments * frag_theo_mass)) {
	AKANTU_DEBUG_ERROR("The mass of small fragments is wrong!");
	return EXIT_FAILURE;
      }
    }

    /// displace fragments
    displaceElements(model, lim, el_size * 2);
    model.updateResidual();
    model.dump();

    lim += el_size;
  }

  /// check centers
  const Array<Real> & fragment_center = model.getFragmentsCenter();

  Real initial_position = -L / 2. + el_size / 2.;

  for (UInt frag = 0; frag < total_nb_fragment; ++frag) {
    Real theoretical_center = initial_position + el_size * frag;

    UInt f_index = 0;
    while (f_index < total_nb_fragment &&
  	   !Math::are_float_equal(fragment_center(f_index, 0), theoretical_center))
      ++f_index;

    if (f_index == total_nb_fragment) {
      AKANTU_DEBUG_ERROR("The fragments' center is wrong!");
      return EXIT_FAILURE;
    }
  }

  // /// Main loop
  // for (UInt s = 1; s <= max_steps; ++s) {

  //   model.checkCohesiveStress();

  //   model.explicitPred();
  //   model.updateResidual();
  //   model.updateAcceleration();
  //   model.explicitCorr();

  //   /// apply boundary conditions
  //   try {
  //     model.applyBC(BC::Dirichlet::IncrementValue(-disp_increment, BC::_x), "Left_side");
  //   } catch(...) {}
  //   try {
  //     model.applyBC(BC::Dirichlet::IncrementValue( disp_increment, BC::_x), "Right_side");
  //   } catch(...) {}

  //   if(s % 1 == 0) {
  //     // model.dump();

  //     std::cout << "passing step " << s << "/" << max_steps << std::endl;

  //     model.computeFragmentsData();

  //     /// check number of fragments
  //     UInt nb_fragment_num = model.getNbFragment();
  //     UInt nb_fragment = countFragments(model);

  //     if (nb_fragment != nb_fragment_num) {
  // 	std::cout << "The number of fragments is wrong at step "
  // 		  << s << std::endl;
  // 	std::cout << "Theoretical: " << nb_fragment
  // 		  << " Computed: " << nb_fragment_num << std::endl;
  // 	return EXIT_FAILURE;
  //     }

  //     /// check mass computation
  //     Real total_mass = 0.;
  //     for (UInt frag = 0; frag < nb_fragment_num; ++frag) {
  // 	total_mass += fragment_mass(frag, 0);
  //     }

  //     if (!Math::are_float_equal(theoretical_mass, total_mass)) {
  // 	std::cout << "The fragments' mass is wrong!" << std::endl;
  // 	return EXIT_FAILURE;
  //     }

  //   }
  // }

  // /// check that all cohesive elements are broken
  // UInt nb_fragment = model.getNbFragment();
  // UInt nb_element = mesh.getNbElement(type);
  // comm.allReduce(&nb_element, 1, _so_sum);

  // if (nb_fragment != nb_element) {
  //   std::cout << "The bar isn't breaking everywhere, increase strain rate" << std::endl;
  //   return EXIT_FAILURE;
  // }

  // /// check velocities
  // const Array<Real> & fragment_velocity = model.getFragmentsVelocity();
  // const Array<Real> & fragment_center = model.getFragmentsCenter();

  // Real fragment_length = L / nb_fragment;
  // Real initial_position = -L / 2. + fragment_length / 2.;

  // for (UInt frag = 0; frag < nb_fragment; ++frag) {
  //   Real theoretical_center = initial_position + fragment_length * frag;

  //   UInt f_index = 0;
  //   while (f_index < nb_fragment &&
  // 	   !Math::are_float_equal(fragment_center(f_index, 0), theoretical_center))
  //     ++f_index;

  //   if (f_index == nb_fragment) {
  //     std::cout << "The fragments' center is wrong!" << std::endl;
  //     return EXIT_FAILURE;
  //   }

  //   Real initial_vel = fragment_center(frag, 0) * strain_rate;

  //   Math::setTolerance(100);

  //   if (!Math::are_float_equal(fragment_velocity(frag), initial_vel)) {
  //     std::cout << "The fragments' velocity is wrong!" << std::endl;
  //     return EXIT_FAILURE;
  //   }
  // }

  finalize();

  if (prank == 0)
    std::cout << "OK: test_cohesive_buildfragments was passed!" << std::endl;
  return EXIT_SUCCESS;
}

void verticalInsertionLimit(SolidMechanicsModelCohesive & model) {
  UInt spatial_dimension = model.getSpatialDimension();
  const Mesh & mesh_facets = model.getMeshFacets();
  const Array<Real> & position = mesh_facets.getNodes();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end();
       ++gt) {
    GhostType ghost_type = *gt;

    Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, ghost_type);
    Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, ghost_type);
    for(; it != end; ++it) {
      ElementType type = *it;

      Array<bool> & check_facets
	= model.getElementInserter().getCheckFacets(type, ghost_type);

      const Array<UInt> & connectivity = mesh_facets.getConnectivity(type, ghost_type);
      UInt nb_nodes_per_facet = connectivity.getNbComponent();

      for (UInt f = 0; f < check_facets.getSize(); ++f) {
	if (!check_facets(f)) continue;

	UInt nb_aligned_nodes = 1;
	Real first_node_pos = position(connectivity(f, 0), 0);

	for (; nb_aligned_nodes < nb_nodes_per_facet; ++nb_aligned_nodes) {
	  Real other_node_pos = position(connectivity(f, nb_aligned_nodes), 0);

	  if (std::abs(first_node_pos - other_node_pos) > 1.e-10)
	    break;
	}

	if (nb_aligned_nodes != nb_nodes_per_facet) {
	  Vector<Real> barycenter(spatial_dimension);

	  mesh_facets.getBarycenter(f, type, barycenter.storage(), ghost_type);
	  std::cout << barycenter(0) << " " << barycenter(1) << std::endl;

	  check_facets(f) = false;
	}
      }
    }
  }
}

void displaceElements(SolidMechanicsModelCohesive & model,
		      const Real lim,
		      const Real amount) {
  UInt spatial_dimension = model.getSpatialDimension();
  Array<Real> & displacement = model.getDisplacement();
  Mesh & mesh = model.getMesh();
  UInt nb_nodes = mesh.getNbNodes();
  Array<bool> displaced(nb_nodes);
  displaced.clear();
  Vector<Real> barycenter(spatial_dimension);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end();
       ++gt) {
    GhostType ghost_type = *gt;

    Mesh::type_iterator it  = mesh.firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, ghost_type);
    for(; it != end; ++it) {
      ElementType type = *it;
      const Array<UInt> & connectivity = mesh.getConnectivity(type, ghost_type);
      UInt nb_element = connectivity.getSize();
      UInt nb_nodes_per_element = connectivity.getNbComponent();

      Array<UInt>::const_iterator<Vector<UInt> > conn_el
	= connectivity.begin(nb_nodes_per_element);

      for (UInt el = 0; el < nb_element; ++el) {
	mesh.getBarycenter(el, type, barycenter.storage(), ghost_type);

	if (barycenter(0) < lim) {
	  const Vector<UInt> & conn = conn_el[el];

	  for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	    UInt node = conn(n);
	    if (!displaced(node)) {
	      displacement(node, 0) -= amount;
	      displaced(node) = true;
	    }
	  }
	}
      }
    }
  }
}
