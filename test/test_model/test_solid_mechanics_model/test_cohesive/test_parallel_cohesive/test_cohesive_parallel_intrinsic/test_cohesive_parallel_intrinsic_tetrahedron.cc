/**
 * @file   test_cohesive_parallel_intrinsic_tetrahedron.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Fri Sep 20 15:59:40 2013
 *
 * @brief  Test for 3D intrinsic cohesive elements simulation in parallel
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
#include "mesh_io.hh"
#include "mesh_utils.hh"
#include "model.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "dumper_paraview.hh"
#include "static_communicator.hh"
#include "dof_synchronizer.hh"
#include "material_cohesive.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

void updateDisplacement(SolidMechanicsModelCohesive & model,
			const ByElementTypeUInt & elements,
			Vector<Real> & increment);

bool checkTractions(SolidMechanicsModelCohesive & model,
		    Vector<Real> & opening,
		    Vector<Real> & theoretical_traction,
		    Matrix<Real> & rotation);

void findNodesToCheck(const Mesh & mesh,
		      const ByElementTypeUInt & elements,
		      Array<UInt> & nodes_to_check);

bool checkEquilibrium(const Mesh & mesh,
		      const Array<Real> & residual);

bool checkResidual(const Array<Real> & residual,
		   const Vector<Real> & traction,
		   const Array<UInt> & nodes_to_check,
		   const Matrix<Real> & rotation);

void findElementsToDisplace(const Mesh & mesh,
			    ByElementTypeUInt & elements);

int main(int argc, char *argv[]) {
  initialize(argc, argv);
  debug::setDebugLevel(dblInfo);

  const UInt spatial_dimension = 3;
  const UInt max_steps = 60;
  const Real increment_constant = 0.01;

  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  UInt nb_nodes_to_check_serial = 0;
  UInt total_nb_nodes = 0;
  Array<UInt> nodes_to_check_serial;

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    // Read the mesh
    mesh.read("tetrahedron.msh");

    /// count nodes with zero position
    const Array<Real> & position = mesh.getNodes();
    for (UInt n = 0; n < position.getSize(); ++n) {
      if (Math::are_float_equal(position(n, 0), 0.))
	++nb_nodes_to_check_serial;
    }

    /// insert cohesive elements
    Array<Real> limits(spatial_dimension, 2);
    limits(0, 0) = -0.01;
    limits(0, 1) = 0.01;
    limits(1, 0) = -100;
    limits(1, 1) = 100;
    limits(2, 0) = -100;
    limits(2, 1) = 100;

    MeshUtils::insertIntrinsicCohesiveElementsInArea(mesh, limits);

    total_nb_nodes = mesh.getNbNodes();

    /// find nodes to check in serial
    ByElementTypeUInt elements_serial("elements_serial", "");
    findElementsToDisplace(mesh, elements_serial);
    findNodesToCheck(mesh, elements_serial, nodes_to_check_serial);

    /// partition the mesh
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    debug::setDebugLevel(dblDump);
    partition->partitionate(psize);
    debug::setDebugLevel(dblInfo);
  }

  comm.broadcast(&nb_nodes_to_check_serial, 1, 0);

  SolidMechanicsModelCohesive model(mesh);

  model.initParallel(partition);
  model.initFull("material_tetrahedron.dat");

  {
    comm.broadcast(&total_nb_nodes, 1, 0);

    Array<Int> nb_local_nodes(psize);
    nb_local_nodes.clear();

    for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
      if (mesh.isLocalOrMasterNode(n)) ++nb_local_nodes(prank);
    }

    comm.allGather(nb_local_nodes.storage(), 1);

    UInt total_nb_nodes_parallel = std::accumulate(nb_local_nodes.begin(),
						   nb_local_nodes.end(), 0);

    Array<UInt> global_nodes_list(total_nb_nodes_parallel);

    UInt first_global_node = std::accumulate(nb_local_nodes.begin(),
					     nb_local_nodes.begin() + prank, 0);

    for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
      if (mesh.isLocalOrMasterNode(n)) {
	global_nodes_list(first_global_node) = mesh.getNodeGlobalId(n);
	++first_global_node;
      }
    }

    comm.allGatherV(global_nodes_list.storage(), nb_local_nodes.storage());

    std::cout << "Maximum node index: "
	      << *(std::max_element(global_nodes_list.begin(),
				    global_nodes_list.end())) << std::endl;

    Array<UInt> repeated_nodes;
    repeated_nodes.resize(0);

    for (UInt n = 0; n < total_nb_nodes_parallel; ++n) {
      UInt appearances = std::count(global_nodes_list.begin() + n,
				    global_nodes_list.end(),
				    global_nodes_list(n));

      if (appearances > 1) {
	std::cout << "Node " << global_nodes_list(n)
		  << " appears " << appearances << " times" << std::endl;

	std::cout << "  in position: " << n;

	repeated_nodes.push_back(global_nodes_list(n));

	UInt * node_position = global_nodes_list.storage() + n;

	for (UInt i = 1; i < appearances; ++i) {
	  node_position = std::find(node_position + 1,
				    global_nodes_list.storage()
				    + total_nb_nodes_parallel,
				    global_nodes_list(n));

	  UInt current_index = node_position - global_nodes_list.storage();

	  std::cout << ", " << current_index;
	}
	std::cout << std::endl << std::endl;
      }
    }


    for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
      UInt global_node = mesh.getNodeGlobalId(n);

      if (std::find(repeated_nodes.begin(), repeated_nodes.end(), global_node)
	  != repeated_nodes.end()) {
	std::cout << "Repeated global node " << global_node
		  << " corresponds to local node " << n << std::endl;
      }
    }


    if (total_nb_nodes != total_nb_nodes_parallel) {
      if (prank == 0) {
	std::cout << "Error: total number of nodes is wrong in parallel" << std::endl;
	std::cout << "Serial: " << total_nb_nodes
		  << " Parallel: " << total_nb_nodes_parallel << std::endl;
      }
      finalize();
      return EXIT_FAILURE;
    }
  }

  model.updateResidual();

  model.setBaseName("intrinsic_parallel_tetrahedron");
  model.addDumpFieldVector("displacement");
  model.addDumpField("residual");
  model.addDumpField("partitions");
  model.dump();

  DumperParaview dumper("cohesive_elements_parallel_tetrahedron");
  dumper.registerMesh(mesh, spatial_dimension, _not_ghost, _ek_cohesive);
  DumperIOHelper::Field * cohesive_displacement =
    new DumperIOHelper::NodalField<Real>(model.getDisplacement());
  cohesive_displacement->setPadding(3);
  dumper.registerField("displacement", cohesive_displacement);

  dumper.dump();

  /// find elements to displace
  ByElementTypeUInt elements("elements", "");
  findElementsToDisplace(mesh, elements);

  /// find nodes to check
  Array<UInt> nodes_to_check;
  findNodesToCheck(mesh, elements, nodes_to_check);

  Vector<Int> nodes_to_check_size(psize);
  nodes_to_check_size(prank) = nodes_to_check.getSize();
  comm.allGather(nodes_to_check_size.storage(), 1);

  UInt nodes_to_check_global_size = std::accumulate(nodes_to_check_size.storage(),
  						    nodes_to_check_size.storage()
  						    + psize, 0);

  Array<UInt> nodes_to_check_global(nodes_to_check_global_size);
  UInt begin_index = std::accumulate(nodes_to_check_size.storage(),
  				     nodes_to_check_size.storage() + prank, 0);

  for (UInt i = begin_index; i < begin_index + nodes_to_check_size(prank); ++i) {
    UInt node = nodes_to_check(i - begin_index);
    nodes_to_check_global(i) = mesh.getNodeGlobalId(node);
  }

  comm.allGatherV(nodes_to_check_global.storage(),
  		  nodes_to_check_size.storage());

  if (nodes_to_check_global_size != nb_nodes_to_check_serial) {
    if (prank == 0) {
      std::cout << "Error: number of nodes to check is wrong in parallel" << std::endl;
      std::cout << "Serial: " << nb_nodes_to_check_serial
  		<< "Parallel: " << nodes_to_check_global_size << std::endl;

      for (UInt n = 0; n < nodes_to_check_serial.getSize(); ++n) {
  	if (std::find(nodes_to_check_global.begin(),
  		      nodes_to_check_global.end(),
  		      nodes_to_check_serial(n)) == nodes_to_check_global.end()) {
  	  std::cout << "Node number " << nodes_to_check_serial(n)
  		    << " not found in parallel" << std::endl;
  	}
      }
    }
    finalize();
    return EXIT_FAILURE;
  }

  /// rotate mesh
  Real angle = 1.;

  Matrix<Real> rotation(spatial_dimension, spatial_dimension);
  rotation.clear();
  rotation(0, 0) = std::cos(angle);
  rotation(0, 1) = std::sin(angle) * -1.;
  rotation(1, 0) = std::sin(angle);
  rotation(1, 1) = std::cos(angle);
  rotation(2, 2) = 1.;


  Vector<Real> increment_tmp(spatial_dimension);
  for (UInt dim = 0; dim < spatial_dimension; ++dim) {
    increment_tmp(dim) = (dim + 1) * increment_constant;
  }

  Vector<Real> increment(spatial_dimension);
  increment.mul<false>(rotation, increment_tmp);

  Array<Real> & position = mesh.getNodes();
  Array<Real> position_tmp(position);

  Array<Real>::iterator<Vector<Real> > position_it = position.begin(spatial_dimension);
  Array<Real>::iterator<Vector<Real> > position_end = position.end(spatial_dimension);
  Array<Real>::iterator<Vector<Real> > position_tmp_it
    = position_tmp.begin(spatial_dimension);

  for (; position_it != position_end; ++position_it, ++position_tmp_it)
    position_it->mul<false>(rotation, *position_tmp_it);

  model.dump();
  dumper.dump();

  updateDisplacement(model, elements, increment);


  Real theoretical_Ed = 0;

  Vector<Real> opening(spatial_dimension);
  Vector<Real> traction(spatial_dimension);
  Vector<Real> opening_old(spatial_dimension);
  Vector<Real> traction_old(spatial_dimension);

  opening.clear();
  traction.clear();
  opening_old.clear();
  traction_old.clear();

  Vector<Real> Dt(spatial_dimension);
  Vector<Real> Do(spatial_dimension);

  const Array<Real> & residual = model.getResidual();

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {

    model.updateResidual();

    opening += increment_tmp;
    if (checkTractions(model, opening, traction, rotation) ||
	checkEquilibrium(mesh, residual) ||
	checkResidual(residual, traction, nodes_to_check, rotation)) {
      finalize();
      return EXIT_FAILURE;
    }

    /// compute energy
    Do = opening;
    Do -= opening_old;

    Dt = traction_old;
    Dt += traction;

    theoretical_Ed += .5 * Do.dot(Dt);

    opening_old = opening;
    traction_old = traction;


    updateDisplacement(model, elements, increment);

    if(s % 10 == 0) {
      std::cout << "passing step " << s << "/" << max_steps << std::endl;
      model.dump();
      dumper.dump();
    }
  }

  model.dump();
  dumper.dump();

  Real Ed = model.getEnergy("dissipated");

  theoretical_Ed *= 4.;

  std::cout << Ed << " " << theoretical_Ed << std::endl;

  if (!Math::are_float_equal(Ed, theoretical_Ed) || std::isnan(Ed)) {

    if (prank == 0)
      std::cout << "Error: the dissipated energy is incorrect" << std::endl;

    finalize();
    return EXIT_FAILURE;
  }

  finalize();
  if(prank == 0) std::cout << "OK: Test passed!" << std::endl;
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */

void updateDisplacement(SolidMechanicsModelCohesive & model,
			const ByElementTypeUInt & elements,
			Vector<Real> & increment) {

  UInt spatial_dimension = model.getSpatialDimension();
  Mesh & mesh = model.getFEM().getMesh();
  UInt nb_nodes = mesh.getNbNodes();

  Array<Real> & displacement = model.getDisplacement();
  Array<bool> update(nb_nodes);
  update.clear();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end();
       ++gt) {

    GhostType ghost_type = *gt;
    Mesh::type_iterator it   = mesh.firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last = mesh.lastType(spatial_dimension, ghost_type);

    for (; it != last; ++it) {
      ElementType type = *it;

      const Array<UInt> & elem = elements(type, ghost_type);
      const Array<UInt> & connectivity = mesh.getConnectivity(type, ghost_type);
      UInt nb_nodes_per_element = connectivity.getNbComponent();

      for (UInt el = 0; el < elem.getSize(); ++el) {
	for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	  UInt node = connectivity(elem(el), n);
	  if (!update(node)) {
	    Vector<Real> node_disp(displacement.storage() + node * spatial_dimension,
				   spatial_dimension);
	    node_disp += increment;
	    update(node) = true;
	  }
	}
      }
    }
  }
}

/* -------------------------------------------------------------------------- */

bool checkTractions(SolidMechanicsModelCohesive & model,
		    Vector<Real> & opening,
		    Vector<Real> & theoretical_traction,
		    Matrix<Real> & rotation) {
  UInt spatial_dimension = model.getSpatialDimension();
  const Mesh & mesh = model.getMesh();

  const MaterialCohesive & mat_cohesive
    = dynamic_cast < const MaterialCohesive & > (model.getMaterial(1));

  const Real sigma_c = mat_cohesive.getParam< RandomInternalField<Real, FacetInternalField> >("sigma_c");
  const Real beta = mat_cohesive.getParam<Real>("beta");
  //  const Real G_cI = mat_cohesive.getParam<Real>("G_cI");
  //  Real G_cII = mat_cohesive.getParam<Real>("G_cII");
  const Real delta_0 = mat_cohesive.getParam<Real>("delta_0");
  const Real kappa = mat_cohesive.getParam<Real>("kappa");
  Real delta_c = delta_0 * sigma_c / (sigma_c - 1.);

  Vector<Real> normal_opening(spatial_dimension);
  normal_opening.clear();
  normal_opening(0) = opening(0);
  Real normal_opening_norm = normal_opening.norm();

  Vector<Real> tangential_opening(spatial_dimension);
  tangential_opening.clear();
  for (UInt dim = 1; dim < spatial_dimension; ++dim)
    tangential_opening(dim) = opening(dim);

  Real tangential_opening_norm = tangential_opening.norm();

  Real beta2_kappa2 = beta * beta/kappa/kappa;
  Real beta2_kappa  = beta * beta/kappa;

  Real delta = std::sqrt(tangential_opening_norm * tangential_opening_norm
			 * beta2_kappa2 +
			 normal_opening_norm * normal_opening_norm);

  delta = std::max(delta, delta_0);

  Real theoretical_damage = std::min(delta / delta_c, 1.);

  if (Math::are_float_equal(theoretical_damage, 1.))
    theoretical_traction.clear();
  else {
    theoretical_traction  = tangential_opening;
    theoretical_traction *= beta2_kappa;
    theoretical_traction += normal_opening;
    theoretical_traction *= sigma_c / delta * (1. - theoretical_damage);
  }

  Vector<Real> theoretical_traction_rotated(spatial_dimension);
  theoretical_traction_rotated.mul<false>(rotation, theoretical_traction);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end();
       ++gt) {

    GhostType ghost_type = *gt;
    Mesh::type_iterator it
      = mesh.firstType(spatial_dimension, ghost_type, _ek_cohesive);
    Mesh::type_iterator last
    = mesh.lastType(spatial_dimension, ghost_type, _ek_cohesive);
    
    for (; it != last; ++it) {
      ElementType type = *it;

      const Array<Real> & traction = mat_cohesive.getTraction(type, ghost_type);
      const Array<Real> & damage = mat_cohesive.getDamage(type, ghost_type);

      UInt nb_quad_per_el
	= model.getFEM("CohesiveFEM").getNbQuadraturePoints(type);
      UInt nb_element = model.getMesh().getNbElement(type, ghost_type);
      UInt tot_nb_quad = nb_element * nb_quad_per_el;

      for (UInt q = 0; q < tot_nb_quad; ++q) {
	for (UInt dim = 0; dim < spatial_dimension; ++dim) {
	  if (!Math::are_float_equal(theoretical_traction_rotated(dim),
				     traction(q, dim))) {
	    std::cout << "Error: tractions are incorrect" << std::endl;
	    return 1;
	  }
	}

	if (ghost_type == _not_ghost)
	  if (!Math::are_float_equal(theoretical_damage, damage(q))) {
	    std::cout << "Error: damage is incorrect" << std::endl;
	    return 1;
	}
      }
    }
  }

  return 0;
}

/* -------------------------------------------------------------------------- */

void findNodesToCheck(const Mesh & mesh,
		      const ByElementTypeUInt & elements,
		      Array<UInt> & nodes_to_check) {
  UInt spatial_dimension = mesh.getSpatialDimension();
  const Array<Real> & position = mesh.getNodes();

  UInt nb_nodes = position.getSize();

  Array<bool> checked_nodes(nb_nodes);
  checked_nodes.clear();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last = mesh.lastType(spatial_dimension);

  for (; it != last; ++it) {
    ElementType type = *it;

    const Array<UInt> & elem = elements(type);
    const Array<UInt> & connectivity = mesh.getConnectivity(type);
    UInt nb_nodes_per_elem = connectivity.getNbComponent();

    for (UInt el = 0; el < elem.getSize(); ++el) {

      UInt element = elem(el);
      Vector<UInt> conn_el(connectivity.storage() + nb_nodes_per_elem * element,
			   nb_nodes_per_elem);

      for (UInt n = 0; n < nb_nodes_per_elem; ++n) {
	UInt node = conn_el(n);
	if (Math::are_float_equal(position(node, 0), 0.)
	    && !checked_nodes(node)
	    && mesh.isLocalOrMasterNode(node)) {
	  checked_nodes(node) = true;
	  nodes_to_check.push_back(node);
	}
      }
    }
  }
}

/* -------------------------------------------------------------------------- */

bool checkEquilibrium(const Mesh & mesh,
		      const Array<Real> & residual) {

  UInt spatial_dimension = residual.getNbComponent();

  Vector<Real> residual_sum(spatial_dimension);
  residual_sum.clear();

  Array<Real>::const_iterator<Vector<Real> > res_it
    = residual.begin(spatial_dimension);

  for (UInt n = 0; n < residual.getSize(); ++n, ++res_it) {
    if (mesh.isLocalOrMasterNode(n))
      residual_sum += *res_it;
  }

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  comm.allReduce(residual_sum.storage(), spatial_dimension, _so_sum);

  for (UInt s = 0; s < spatial_dimension; ++s) {
    if (!Math::are_float_equal(residual_sum(s), 0.)) {

      if (comm.whoAmI() == 0)
	std::cout << "Error: system is not in equilibrium!" << std::endl;

      return 1;
    }
  }

  return 0;
}

/* -------------------------------------------------------------------------- */

bool checkResidual(const Array<Real> & residual,
		   const Vector<Real> & traction,
		   const Array<UInt> & nodes_to_check,
		   const Matrix<Real> & rotation) {

  UInt spatial_dimension = residual.getNbComponent();

  Vector<Real> total_force(spatial_dimension);
  total_force.clear();

  for (UInt n = 0; n < nodes_to_check.getSize(); ++n) {
    UInt node = nodes_to_check(n);

    Vector<Real> res(residual.storage() + node * spatial_dimension,
		     spatial_dimension);

    total_force += res;
  }

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  comm.allReduce(total_force.storage(), spatial_dimension, _so_sum);

  Vector<Real> theoretical_total_force(spatial_dimension);
  theoretical_total_force.mul<false>(rotation, traction);
  theoretical_total_force *= -1 * 2 * 2;

  for (UInt s = 0; s < spatial_dimension; ++s) {
    if (!Math::are_float_equal(total_force(s), theoretical_total_force(s))) {

      if (comm.whoAmI() == 0)
	std::cout << "Error: total force isn't correct!" << std::endl;

      return 1;
    }
  }

  return 0;
}

/* -------------------------------------------------------------------------- */

void findElementsToDisplace(const Mesh & mesh,
			    ByElementTypeUInt & elements) {
  UInt spatial_dimension = mesh.getSpatialDimension();

  mesh.initByElementTypeArray(elements, 1, spatial_dimension);

  Vector<Real> bary(spatial_dimension);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end();
       ++gt) {

    GhostType ghost_type = *gt;
    Mesh::type_iterator it   = mesh.firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last = mesh.lastType(spatial_dimension, ghost_type);

    for (; it != last; ++it) {
      ElementType type = *it;

      Array<UInt> & elem = elements(type, ghost_type);
      UInt nb_element = mesh.getNbElement(type, ghost_type);

      for (UInt el = 0; el < nb_element; ++el) {
	mesh.getBarycenter(el, type, bary.storage(), ghost_type);
	if (bary(0) > 0.01) elem.push_back(el);
      }
    }
  }
}
