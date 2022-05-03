/**
 * @file   RVE_tools.cc
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Tue Jan 16 10:26:53 2014
 * @update Tue Feb 8  2022
 *
 * @brief  implementation tools for the analysis of RVEs
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
 **/

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/property_map/property_map.hpp>
#include <fstream>
/* -------------------------------------------------------------------------- */
#include "aka_voigthelper.hh"
#include "communicator.hh"
#include "crack_numbers_updater.hh"
#include "material_FE2.hh"
#include "material_cohesive_linear_sequential.hh"
#include "material_damage_iterative_orthotropic.hh"
#include "material_iterative_stiffness_reduction.hh"
#include "mesh_utils.hh"
#include "node_synchronizer.hh"
#include "nodes_eff_stress_updater.hh"
#include "nodes_flag_updater.hh"
#include "non_linear_solver.hh"
#include "rve_tools.hh"
#include "solid_mechanics_model.hh"
#include "solid_mechanics_model_RVE.hh"
#include "solid_mechanics_model_cohesive.hh"
#include <cmath>
#include <mesh_events.hh>

/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
RVETools::RVETools(SolidMechanicsModel & model)
    : model(model), volume(0.), nb_dumps(0), cohesive_insertion(false),
      doubled_facets_ready(false), crack_central_nodes_ready(false),
      disp_stored(0, model.getSpatialDimension()),
      ext_force_stored(0, model.getSpatialDimension()),
      boun_stored(0, model.getSpatialDimension()),
      tensile_homogenization(false), modified_pos(false),
      partition_border_nodes(false), nodes_eff_stress(0), crack_nodes(false) {

  // register event handler for rve tools
  auto & mesh = model.getMesh();
  mesh.registerEventHandler(*this, _ehp_synchronizer);
  // mesh.registerEventHandler(*this, _ehp_lowest);
  if (mesh.hasMeshFacets()) {
    mesh.getMeshFacets().registerEventHandler(*this, _ehp_synchronizer);
    // mesh.getMeshFacets().registerEventHandler(*this, _ehp_lowest);
  }

  // find four corner nodes of the RVE
  findCornerNodes();

  // initialize the modified_pos array with the initial number of nodes
  auto nb_nodes = mesh.getNbNodes();
  modified_pos.resize(nb_nodes);
  modified_pos.set(false);
  partition_border_nodes.resize(nb_nodes);
  partition_border_nodes.set(false);
  nodes_eff_stress.resize(nb_nodes);
  nodes_eff_stress.set(0);
  crack_nodes.resize(nb_nodes);
  crack_nodes.set(false);
}
/* -------------------------------------------------------------------------- */
void RVETools::computePhaseVolumes() {
  // compute volume of each phase and save it into a map
  for (auto && mat : model.getMaterials()) {
    this->phase_volumes[mat.getName()] = computePhaseVolume(mat.getName());
    auto it = this->phase_volumes.find(mat.getName());
    if (it == this->phase_volumes.end()) {
      this->phase_volumes.erase(mat.getName());
    }
  }
}
/* -------------------------------------------------------------------------- */
void RVETools::computeModelVolume() {
  auto const dim = model.getSpatialDimension();
  auto & mesh = model.getMesh();
  auto & fem = model.getFEEngine("SolidMechanicsFEEngine");
  GhostType gt = _not_ghost;
  this->volume = 0;
  for (auto element_type : mesh.elementTypes(dim, gt, _ek_regular)) {
    Array<Real> Volume(mesh.getNbElement(element_type) *
                           fem.getNbIntegrationPoints(element_type),
                       1, 1.);
    this->volume += fem.integrate(Volume, element_type);
  }

  // do not communicate if model is multi-scale
  if (not aka::is_of_type<SolidMechanicsModelRVE>(model)) {
    auto && comm = akantu::Communicator::getWorldCommunicator();
    comm.allReduce(this->volume, SynchronizerOperation::_sum);
  }
}
/* ------------------------------------------------------------------------- */
void RVETools::applyFreeExpansionBC(Real offset) {
  // boundary conditions
  const auto & mesh = model.getMesh();
  const Vector<Real> & lowerBounds = mesh.getLowerBounds();
  const Vector<Real> & upperBounds = mesh.getUpperBounds();
  const UInt dim = mesh.getSpatialDimension();
  const Array<Real> & pos = mesh.getNodes();
  Array<Real> & disp = model.getDisplacement();
  Array<bool> & boun = model.getBlockedDOFs();

  // accessing bounds
  Real left = lowerBounds(0) + offset;
  Real right = upperBounds(0) - offset;

  switch (dim) {
  case 2: {
    Real bottom = lowerBounds(1) + offset;
    Real top = upperBounds(1) - offset;
    Real eps = std::abs((top - bottom) * 1e-6);
    for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
      if ((left <= pos(i, 0)) and (pos(i, 0) <= right)) {
        if (std::abs(pos(i, 1) - bottom) < eps) {
          boun(i, 1) = true;
          disp(i, 1) = 0.0;

          if (std::abs(pos(i, 0) - left) < eps) {
            boun(i, 0) = true;
            disp(i, 0) = 0.0;
          }
        }
      }
    }
    break;
  }
  case 3: {
    // accessing bounds
    Real bottom = lowerBounds(2) + offset;
    Real top = upperBounds(2) - offset;
    Real back = lowerBounds(1) + offset;
    Real front = upperBounds(1) - offset;
    Real eps = std::abs((top - bottom) * 1e-6);
    for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
      if ((left <= pos(i, 0)) and (pos(i, 0) <= right) and
          (bottom <= pos(i, 2)) and (pos(i, 2) <= top) and
          (back <= pos(i, 1)) and (pos(i, 1) <= front)) {

        if (std::abs(pos(i, 2) - bottom) < eps) {
          boun(i, 2) = true;
          disp(i, 2) = 0.0;
        }
        if (std::abs(pos(i, 0) - left) < eps) {
          boun(i, 0) = true;
          disp(i, 0) = 0.0;
        }
        if (std::abs(pos(i, 1) - back) < eps) {
          boun(i, 1) = true;
          disp(i, 1) = 0.0;
        }
      }
    }
    break;
  }
  }
}

/* ------------------------------------------------------------------- */
void RVETools::applyLoadedBC(const Vector<Real> & traction,
                             const ID & element_group, bool multi_axial) {

  // boundary conditions
  const auto & mesh = model.getMesh();
  const UInt dim = mesh.getSpatialDimension();
  const Vector<Real> & lowerBounds = mesh.getLowerBounds();
  const Vector<Real> & upperBounds = mesh.getUpperBounds();
  Real bottom = lowerBounds(1);
  Real left = lowerBounds(0);
  Real top = upperBounds(1);

  Real eps = std::abs((top - bottom) * 1e-6);
  const Array<Real> & pos = mesh.getNodes();
  Array<Real> & disp = model.getDisplacement();
  Array<bool> & boun = model.getBlockedDOFs();

  switch (dim) {
  case 2: {
    for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
      if (std::abs(pos(i, 1) - bottom) < eps) {
        boun(i, 1) = true;
        disp(i, 1) = 0.0;
      }
      if ((std::abs(pos(i, 1) - bottom) < eps) &&
          (std::abs(pos(i, 0) - left) < eps)) {
        boun(i, 0) = true;
        disp(i, 0) = 0.0;
      }
      if (multi_axial && (std::abs(pos(i, 0) - left) < eps)) {
        boun(i, 0) = true;
        disp(i, 0) = 0.0;
      }
    }
    break;
  }
  case 3: {
    Real back = lowerBounds(2);
    for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
      if ((std::abs(pos(i, 1) - bottom) < eps) &&
          (std::abs(pos(i, 0) - left) < eps) &&
          (std::abs(pos(i, 2) - back) < eps)) {
        boun(i, 0) = true;
        boun(i, 1) = true;
        boun(i, 2) = true;
        disp(i, 0) = 0.0;
        disp(i, 1) = 0.0;
        disp(i, 2) = 0.0;
      }
      if (std::abs(pos(i, 1) - bottom) < eps) {
        boun(i, 1) = true;
        disp(i, 1) = 0.0;
      }
      if (multi_axial && (std::abs(pos(i, 2) - back) < eps)) {
        boun(i, 2) = true;
        disp(i, 2) = 0.0;
      }
      if (multi_axial && (std::abs(pos(i, 0) - left) < eps)) {
        boun(i, 0) = true;
        disp(i, 0) = 0.0;
      }
    }
  }
  }
  // try {
  model.applyBC(BC::Neumann::FromTraction(traction), element_group);
  // } catch (...) {
  // }
}

/* ------------------------------------------------------------------- */
void RVETools::applyExternalTraction(Real traction,
                                     SpatialDirection direction) {

  // boundary conditions
  auto && mesh = model.getMesh();
  auto && upperBounds = mesh.getUpperBounds();
  auto limit = upperBounds(direction);

  auto && pos = mesh.getNodes();
  auto & force = model.getExternalForce();

  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    if (Math::are_float_equal(pos(i, direction), limit)) {
      force(i, direction) = traction;
    }
  }
}

/* --------------------------------------------------------------------------
 */
void RVETools::fillNodeGroup(NodeGroup & node_group, SpatialDirection dir,
                             Real offset) {
  const auto & mesh = model.getMesh();
  const Vector<Real> & upperBounds = mesh.getUpperBounds();
  const Vector<Real> & lowerBounds = mesh.getLowerBounds();
  Real top = upperBounds(dir) - offset;
  Real bottom = lowerBounds(dir) + offset;
  Real eps = std::abs((top - bottom) * 1e-6);
  const Array<Real> & pos = mesh.getNodes();
  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    if (std::abs(pos(i, dir) - top) < eps) {
      node_group.add(i);
    }
  }
}

/* ------------------------------------------------------------------- */
Real RVETools::computeAverageStrain(SpatialDirection direction, Real offset) {

  Real av_displ{0};
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  const Vector<Real> & upperBounds = mesh.getUpperBounds();
  const Vector<Real> & lowerBounds = mesh.getLowerBounds();
  Real right = upperBounds(0) - offset;
  Real left = lowerBounds(0) + offset;
  Real top = upperBounds(1) - offset;
  Real bottom = lowerBounds(1) + offset;
  Real front;
  Real back;
  Real length{0};
  if (dim == 3) {
    front = upperBounds(2) - offset;
    back = lowerBounds(2) + offset;
  }
  Real eps = std::abs((right - left) * 1e-6);
  const Array<Real> & pos = mesh.getNodes();
  const Array<Real> & disp = model.getDisplacement();
  UInt nb_nodes = 0;

  if (direction == _x) {
    length = std::abs(right - left);
    for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
      if ((bottom <= pos(i, 1)) and (pos(i, 1) <= top) and
          (std::abs(pos(i, 0) - right) < eps) and mesh.isLocalOrMasterNode(i)) {
        if ((dim == 3) and (back <= pos(i, 2)) and (pos(i, 2) <= front)) {
          // dirtly trick to exclude blowed up elements
          if (disp(i, 0) < length / 50) {
            av_displ += disp(i, 0);
            ++nb_nodes;
          }
        }
      }
    }
  }

  else if (direction == _y) {
    length = std::abs(top - bottom);
    for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
      if ((left <= pos(i, 0)) and (pos(i, 0) <= right) and
          (std::abs(pos(i, 1) - top) < eps) and mesh.isLocalOrMasterNode(i)) {
        if ((dim == 3) and (back <= pos(i, 2)) and (pos(i, 2) <= front)) {
          if (disp(i, 1) < length / 50) {
            av_displ += disp(i, 1);
            ++nb_nodes;
          }
        }
      }
    }
  } else if ((direction == _z) && (model.getSpatialDimension() == 3)) {
    AKANTU_DEBUG_ASSERT(model.getSpatialDimension() == 3,
                        "no z-direction in 2D problem");
    length = std::abs(front - back);
    for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
      if ((left <= pos(i, 0)) and (pos(i, 0) <= right) and
          (bottom <= pos(i, 1)) and (pos(i, 1) <= top) and
          (std::abs(pos(i, 2) - front) < eps) and mesh.isLocalOrMasterNode(i)) {
        if (disp(i, 2) < length / 50) {
          av_displ += disp(i, 2);
          ++nb_nodes;
        }
      }
    }
  } else {
    AKANTU_EXCEPTION("The parameter for the testing direction is wrong!!!");
  }

  // do not communicate if model is multi-scale
  if (!aka::is_of_type<SolidMechanicsModelRVE>(model)) {
    auto && comm = akantu::Communicator::getWorldCommunicator();
    comm.allReduce(av_displ, SynchronizerOperation::_sum);
    comm.allReduce(nb_nodes, SynchronizerOperation::_sum);
  }

  return av_displ / nb_nodes / length;
}

/* --------------------------------------------------------------------------
 */
Real RVETools::computeVolumetricExpansion(SpatialDirection direction) {
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  UInt tensor_component = 0;
  if (dim == 2) {
    if (direction == _x) {
      tensor_component = 0;
    } else if (direction == _y) {
      tensor_component = 3;
    } else {
      AKANTU_EXCEPTION("The parameter for the testing direction is wrong!!!");
    }
  }

  else if (dim == 3) {
    if (direction == _x) {
      tensor_component = 0;
    } else if (direction == _y) {
      tensor_component = 4;
    } else if (direction == _z) {
      tensor_component = 8;
    } else {
      AKANTU_EXCEPTION("The parameter for the testing direction is wrong!!!");
    }
  } else {
    AKANTU_EXCEPTION("This example does not work for 1D!!!");
  }

  Real gradu_tot = 0;
  Real tot_volume = 0;
  GhostType gt = _not_ghost;

  for (auto element_type : mesh.elementTypes(dim, gt, _ek_regular)) {
    const FEEngine & fe_engine = model.getFEEngine();
    for (UInt m = 0; m < model.getNbMaterials(); ++m) {
      const ElementTypeMapUInt & element_filter_map =
          model.getMaterial(m).getElementFilter();
      if (!element_filter_map.exists(element_type, gt)) {
        continue;
      }
      const Array<UInt> & elem_filter =
          model.getMaterial(m).getElementFilter(element_type);
      if (elem_filter.empty()) {
        continue;
      }
      Array<Real> gradu_vec(elem_filter.size() *
                                fe_engine.getNbIntegrationPoints(element_type),
                            dim * dim, 0.);

      fe_engine.gradientOnIntegrationPoints(model.getDisplacement(), gradu_vec,
                                            dim, element_type, gt, elem_filter);

      Array<Real> int_gradu_vec(elem_filter.size(), dim * dim, "int_of_gradu");

      fe_engine.integrate(gradu_vec, int_gradu_vec, dim * dim, element_type,
                          _not_ghost, elem_filter);

      for (UInt k = 0; k < elem_filter.size(); ++k) {
        gradu_tot += int_gradu_vec(k, tensor_component);
      }
    }

    Array<Real> Volume(mesh.getNbElement(element_type) *
                           fe_engine.getNbIntegrationPoints(element_type),
                       1, 1.);
    Real int_volume = fe_engine.integrate(Volume, element_type);
    tot_volume += int_volume;
  }

  // do not communicate if model is multi-scale
  if (not aka::is_of_type<SolidMechanicsModelRVE>(model)) {
    auto && comm = akantu::Communicator::getWorldCommunicator();
    comm.allReduce(gradu_tot, SynchronizerOperation::_sum);
    comm.allReduce(tot_volume, SynchronizerOperation::_sum);
  }

  return gradu_tot / tot_volume;
}

/* --------------------------------------------------------------------------
 */
Real RVETools::computeDamagedVolume(const ID & mat_name) {
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  Real total_damage = 0;
  GhostType gt = _not_ghost;
  const Material & mat = model.getMaterial(mat_name);

  for (auto element_type : mesh.elementTypes(dim, gt, _ek_not_defined)) {
    const FEEngine & fe_engine = model.getFEEngine();
    const ElementTypeMapUInt & element_filter_map = mat.getElementFilter();
    if (!element_filter_map.exists(element_type, gt)) {
      continue;
    }
    const Array<UInt> & elem_filter = mat.getElementFilter(element_type);
    if (elem_filter.empty()) {
      continue;
    }
    const Array<Real> & damage = mat.getInternal<Real>("damage")(element_type);

    total_damage += fe_engine.integrate(damage, element_type, gt, elem_filter);
  }

  // do not communicate if model is multi-scale
  if (not aka::is_of_type<SolidMechanicsModelRVE>(model)) {
    auto && comm = akantu::Communicator::getWorldCommunicator();
    comm.allReduce(total_damage, SynchronizerOperation::_sum);
  }

  return total_damage;
}

/* --------------------------------------------------------------------------
 */
void RVETools::computeStiffnessReduction(std::ofstream & file_output, Real time,
                                         bool tension) {

  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  AKANTU_DEBUG_ASSERT(dim != 1, "Example does not work for 1D");

  // save nodal values before test
  storeNodalFields();

  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();

  if (dim == 2) {
    Real int_residual_x = 0;
    Real int_residual_y = 0;
    int_residual_x = performLoadingTest(_x, tension);
    int_residual_y = performLoadingTest(_y, tension);

    if (prank == 0) {
      file_output << time << "," << int_residual_x << "," << int_residual_y
                  << std::endl;
    }
  } else {
    Real int_residual_x = performLoadingTest(_x, tension);
    Real int_residual_y = performLoadingTest(_y, tension);
    Real int_residual_z = performLoadingTest(_z, tension);
    if (prank == 0) {
      file_output << time << "," << int_residual_x << "," << int_residual_y
                  << "," << int_residual_z << std::endl;
    }
  }

  // return the nodal values
  restoreNodalFields();
}

/* -------------------------------------------------------------------------- */
void RVETools::storeNodalFields() {
  auto & disp = this->model.getDisplacement();
  auto & boun = this->model.getBlockedDOFs();
  auto & ext_force = this->model.getExternalForce();
  this->disp_stored.copy(disp);
  this->boun_stored.copy(boun);
  this->ext_force_stored.copy(ext_force);
}

/* -------------------------------------------------------------------------- */
void RVETools::restoreNodalFields() {
  auto & disp = this->model.getDisplacement();
  auto & boun = this->model.getBlockedDOFs();
  auto & ext_force = this->model.getExternalForce();
  disp.copy(this->disp_stored);
  boun.copy(this->boun_stored);
  ext_force.copy(this->ext_force_stored);
}

/* ---------------------------------------------------------------------- */
void RVETools::resetNodalFields() {
  auto & disp = this->model.getDisplacement();
  auto & boun = this->model.getBlockedDOFs();
  auto & ext_force = this->model.getExternalForce();
  disp.zero();
  boun.zero();
  ext_force.zero();
}

/* -------------------------------------------------------------------------- */
void RVETools::restoreInternalFields() {
  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    Material & mat = model.getMaterial(m);
    mat.restorePreviousState();
  }
}
/* -------------------------------------------------------------------------- */
void RVETools::resetInternalFields() {
  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    Material & mat = model.getMaterial(m);
    mat.resetInternalsWithHistory();
  }
}
/* -------------------------------------------------------------------------- */
void RVETools::storeDamageField() {
  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    Material & mat = model.getMaterial(m);
    if (mat.isInternal<Real>("damage_stored", _ek_regular) &&
        mat.isInternal<UInt>("reduction_step_stored", _ek_regular)) {
      for (const auto & el_type : mat.getElementFilter().elementTypes(
               _all_dimensions, _not_ghost, _ek_not_defined)) {
        auto & red_stored =
            mat.getInternal<UInt>("reduction_step_stored")(el_type, _not_ghost);
        auto & red_current =
            mat.getInternal<UInt>("reduction_step")(el_type, _not_ghost);
        red_stored.copy(red_current);

        auto & dam_stored =
            mat.getInternal<Real>("damage_stored")(el_type, _not_ghost);
        auto & dam_current =
            mat.getInternal<Real>("damage")(el_type, _not_ghost);
        dam_stored.copy(dam_current);
      }
    }
  }
}
/* -------------------------------------------------------------------------- */
void RVETools::restoreDamageField() {
  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    Material & mat = model.getMaterial(m);
    if (mat.isInternal<Real>("damage_stored", _ek_regular) &&
        mat.isInternal<UInt>("reduction_step_stored", _ek_regular)) {
      for (const auto & el_type : mat.getElementFilter().elementTypes(
               _all_dimensions, _not_ghost, _ek_not_defined)) {
        auto & red_stored =
            mat.getInternal<UInt>("reduction_step_stored")(el_type, _not_ghost);
        auto & red_current =
            mat.getInternal<UInt>("reduction_step")(el_type, _not_ghost);
        red_current.copy(red_stored);

        auto & dam_stored =
            mat.getInternal<Real>("damage_stored")(el_type, _not_ghost);
        auto & dam_current =
            mat.getInternal<Real>("damage")(el_type, _not_ghost);
        dam_current.copy(dam_stored);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
Real RVETools::performLoadingTest(SpatialDirection direction, bool tension) {
  UInt dir;

  if (direction == _x) {
    dir = 0;
  }
  if (direction == _y) {
    dir = 1;
  }
  if (direction == _z) {
    AKANTU_DEBUG_ASSERT(model.getSpatialDimension() == 3,
                        "Error in problem dimension!!!");
    dir = 2;
  }

  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();

  // indicate to material that its in the loading test
  UInt nb_materials = model.getNbMaterials();
  for (UInt m = 0; m < nb_materials; ++m) {
    Material & mat = model.getMaterial(m);
    if (aka::is_of_type<MaterialDamageIterativeOrthotropic<2>>(mat)) {
      mat.setParam("loading_test", true);
    }
  }

  // boundary conditions
  const Vector<Real> & lowerBounds = mesh.getLowerBounds();
  const Vector<Real> & upperBounds = mesh.getUpperBounds();
  Real bottom = lowerBounds(1);
  Real top = upperBounds(1);
  Real left = lowerBounds(0);
  Real back;

  if (dim == 3) {
    back = lowerBounds(2);
  }

  Real eps = std::abs((top - bottom) * 1e-6);
  Real imposed_displacement = std::abs(lowerBounds(dir) - upperBounds(dir));
  imposed_displacement *= 0.001;
  const auto & pos = mesh.getNodes();
  auto & disp = model.getDisplacement();
  auto & boun = model.getBlockedDOFs();
  auto & ext_force = model.getExternalForce();
  UInt nb_nodes = mesh.getNbNodes();

  disp.zero();
  boun.zero();
  ext_force.zero();

  if (dim == 3) {
    for (UInt i = 0; i < nb_nodes; ++i) {

      // fix one corner node to avoid sliding
      if ((std::abs(pos(i, 1) - bottom) < eps) &&
          (std::abs(pos(i, 0) - left) < eps) &&
          (std::abs(pos(i, 2) - back) < eps)) {
        boun(i, 0) = true;
        boun(i, 1) = true;
        boun(i, 2) = true;
        disp(i, 0) = 0.0;
        disp(i, 1) = 0.0;
        disp(i, 2) = 0.0;
      }

      if ((std::abs(pos(i, dir) - lowerBounds(dir)) < eps)) {
        boun(i, dir) = true;
        disp(i, dir) = 0.0;
      }
      if ((std::abs(pos(i, dir) - upperBounds(dir)) < eps)) {
        boun(i, dir) = true;
        disp(i, dir) = (2 * tension - 1) * imposed_displacement;
      }
    }
  } else {
    for (UInt i = 0; i < nb_nodes; ++i) {

      // fix one corner node to avoid sliding
      if ((std::abs(pos(i, 1) - bottom) < eps) &&
          (std::abs(pos(i, 0) - left) < eps)) {
        boun(i, 0) = true;
        boun(i, 1) = true;
        disp(i, 0) = 0.0;
        disp(i, 1) = 0.0;
      }

      if ((std::abs(pos(i, dir) - lowerBounds(dir)) < eps)) {
        boun(i, dir) = true;
        disp(i, dir) = 0.0;
      }
      if ((std::abs(pos(i, dir) - upperBounds(dir)) < eps)) {
        boun(i, dir) = true;
        disp(i, dir) = (2 * tension - 1) * imposed_displacement;
      }
    }
  }

  try {
    model.solveStep();
  } catch (debug::Exception & e) {
    auto & solver = model.getNonLinearSolver("static");
    int nb_iter = solver.get("nb_iterations");
    std::cout << "Loading test did not converge in " << nb_iter
              << " iterations." << std::endl;
    throw e;
  }

  // compute the force (residual in this case) along the edge of the
  // imposed displacement
  Real int_residual = 0.;
  const Array<Real> & residual = model.getInternalForce();

  for (UInt n = 0; n < nb_nodes; ++n) {
    if (std::abs(pos(n, dir) - upperBounds(dir)) < eps &&
        mesh.isLocalOrMasterNode(n)) {
      int_residual += -residual(n, dir);
    }
  }

  // indicate to material that its out of the loading test
  for (UInt m = 0; m < nb_materials; ++m) {
    Material & mat = model.getMaterial(m);
    if (aka::is_of_type<MaterialDamageIterativeOrthotropic<2>>(mat)) {
      mat.setParam("loading_test", false);
    }
  }

  // restore historical internal fields
  restoreInternalFields();

  // do not communicate if model is multi-scale
  if (not aka::is_of_type<SolidMechanicsModelRVE>(model)) {
    auto && comm = akantu::Communicator::getWorldCommunicator();
    comm.allReduce(int_residual, SynchronizerOperation::_sum);
  }

  return int_residual;
}

/* --------------------------------------------------------------------------
 */
void RVETools::computeAverageProperties(std::ofstream & file_output) {

  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  AKANTU_DEBUG_ASSERT(dim != 1, "Example does not work for 1D");

  Real av_tens_strain_x = computeVolumetricExpansion(_x);
  Real av_tens_strain_y = computeVolumetricExpansion(_y);
  Real av_lin_strain_x = computeAverageStrain(_x);
  Real av_lin_strain_y = computeAverageStrain(_y);
  Real damage_agg, damage_paste, crack_agg, crack_paste;
  computeDamageRatioPerMaterial(damage_agg, "aggregate");
  computeDamageRatioPerMaterial(damage_paste, "paste");
  computeCrackVolumePerMaterial(crack_agg, "aggregate");
  computeCrackVolumePerMaterial(crack_paste, "paste");

  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();

  if (dim == 2) {

    if (prank == 0)
      file_output << av_tens_strain_x << "," << av_tens_strain_y << ","
                  << av_lin_strain_x << "," << av_lin_strain_y << ","
                  << damage_agg << "," << damage_paste << "," << crack_agg
                  << "," << crack_paste << std::endl;
  }

  else {
    Real av_lin_strain_z = computeAverageStrain(_z);
    Real av_tens_strain_z = computeVolumetricExpansion(_z);

    if (prank == 0)
      file_output << av_tens_strain_x << "," << av_tens_strain_y << ","
                  << av_tens_strain_z << "," << av_lin_strain_x << ","
                  << av_lin_strain_y << "," << av_lin_strain_z << ","
                  << damage_agg << "," << damage_paste << "," << crack_agg
                  << "," << crack_paste << std::endl;
  }
}

/* --------------------------------------------------------------------------
 */
void RVETools::computeAverageProperties(std::ofstream & file_output,
                                        Real time) {
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();

  AKANTU_DEBUG_ASSERT(dim != 1, "Example does not work for 1D");

  Real av_tens_strain_x = computeVolumetricExpansion(_x);
  Real av_tens_strain_y = computeVolumetricExpansion(_y);
  Real av_lin_strain_x = computeAverageStrain(_x);
  Real av_lin_strain_y = computeAverageStrain(_y);
  Real damage_agg, damage_paste, crack_agg, crack_paste;
  computeDamageRatioPerMaterial(damage_agg, "aggregate");
  computeDamageRatioPerMaterial(damage_paste, "paste");
  computeCrackVolumePerMaterial(crack_agg, "aggregate");
  computeCrackVolumePerMaterial(crack_paste, "paste");

  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();

  if (dim == 2) {

    if (prank == 0)
      file_output << time << "," << av_tens_strain_x << "," << av_tens_strain_y
                  << "," << av_lin_strain_x << "," << av_lin_strain_y << ","
                  << damage_agg << "," << damage_paste << "," << crack_agg
                  << "," << crack_paste << std::endl;
  }

  else {
    Real av_lin_strain_z = computeAverageStrain(_z);
    Real av_tens_strain_z = computeVolumetricExpansion(_z);

    if (prank == 0)
      file_output << time << "," << av_tens_strain_x << "," << av_tens_strain_y
                  << "," << av_tens_strain_z << "," << av_lin_strain_x << ","
                  << av_lin_strain_y << "," << av_lin_strain_z << ","
                  << damage_agg << "," << damage_paste << "," << crack_agg
                  << "," << crack_paste << std::endl;
  }
}
/* --------------------------------------------------------------------------
 */
void RVETools::computeAveragePropertiesFe2Material(std::ofstream & file_output,
                                                   Real time) {
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();

  AKANTU_DEBUG_ASSERT(dim != 1, "Example does not work for 1D");

  Real av_tens_strain_x = computeVolumetricExpansion(_x);
  Real av_tens_strain_y = computeVolumetricExpansion(_y);
  Real av_lin_strain_x = computeAverageStrain(_x);
  Real av_lin_strain_y = computeAverageStrain(_y);
  Real crack_agg = averageScalarField("crack_volume_ratio_agg");
  Real crack_paste = averageScalarField("crack_volume_ratio_paste");

  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();

  if (dim == 2) {

    if (prank == 0)
      file_output << time << "," << av_tens_strain_x << "," << av_tens_strain_y
                  << "," << av_lin_strain_x << "," << av_lin_strain_y << ","
                  << crack_agg << "," << crack_paste << std::endl;
  }

  else {
    Real av_lin_strain_z = computeAverageStrain(_z);
    Real av_tens_strain_z = computeVolumetricExpansion(_z);

    if (prank == 0)
      file_output << time << "," << av_tens_strain_x << "," << av_tens_strain_y
                  << "," << av_tens_strain_z << "," << av_lin_strain_x << ","
                  << av_lin_strain_y << "," << av_lin_strain_z << ","
                  << crack_agg << "," << crack_paste << std::endl;
  }
}

/* ------------------------------------------------------------------- */
void RVETools::saveState(UInt current_step) {

  // variables for parallel execution
  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();
  GhostType gt = _not_ghost;

  // save the reduction step on each quad
  UInt nb_materials = model.getNbMaterials();

  // open the output file
  std::stringstream file_stream;
  file_stream << "./restart/reduction_steps_file_" << prank << "_"
              << current_step << ".txt";
  std::string reduction_steps_file = file_stream.str();
  std::ofstream file_output;
  file_output.open(reduction_steps_file.c_str());

  if (!file_output.is_open())
    AKANTU_EXCEPTION("Could not create the file " + reduction_steps_file +
                     ", does its folder exist?");

  for (UInt m = 0; m < nb_materials; ++m) {
    const Material & mat = model.getMaterial(m);
    if (mat.getName() == "gel")
      continue;
    // get the reduction steps internal field
    const InternalField<UInt> & reduction_steps =
        mat.getInternal<UInt>("damage_step");
    // loop over all types in that material
    for (auto element_type : reduction_steps.elementTypes(gt)) {
      const Array<UInt> & elem_filter = mat.getElementFilter(element_type, gt);
      if (!elem_filter.size())
        continue;
      const Array<UInt> & reduction_step_array =
          reduction_steps(element_type, gt);
      for (UInt i = 0; i < reduction_step_array.size(); ++i) {
        file_output << reduction_step_array(i) << std::endl;

        if (file_output.fail())
          AKANTU_EXCEPTION("Error in writing data to file " +
                           reduction_steps_file);
      }
    }
  }

  // close the file
  file_output.close();
  comm.barrier();
  // write the number of the last successfully saved step in a file
  if (prank == 0) {
    std::string current_step_file = "./restart/current_step.txt";
    std::ofstream file_output;
    file_output.open(current_step_file.c_str());
    file_output << current_step << std::endl;
    if (file_output.fail())
      AKANTU_EXCEPTION("Error in writing data to file " + current_step_file);
    file_output.close();
  }
}

/* ------------------------------------------------------------------ */
bool RVETools::loadState(UInt & current_step) {
  current_step = 0;
  bool restart = false;

  // variables for parallel execution
  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();
  GhostType gt = _not_ghost;

  // proc 0 has to save the current step
  if (prank == 0) {
    std::string line;
    std::string current_step_file = "./restart/current_step.txt";
    std::ifstream file_input;
    file_input.open(current_step_file.c_str());
    if (file_input.is_open()) {
      std::getline(file_input, line);
      std::stringstream sstr(line);
      sstr >> current_step;
      file_input.close();
    }
  }

  comm.broadcast(current_step);

  if (!current_step)
    return restart;

  if (prank == 0)
    std::cout << "....Restarting simulation" << std::endl;
  restart = true;

  // save the reduction step on each quad
  UInt nb_materials = model.getNbMaterials();

  // open the output file
  std::stringstream file_stream;
  file_stream << "./restart/reduction_steps_file_" << prank << "_"
              << current_step << ".txt";
  std::string reduction_steps_file = file_stream.str();
  std::ifstream file_input;
  file_input.open(reduction_steps_file.c_str());

  if (!file_input.is_open())
    AKANTU_EXCEPTION("Could not open file " + reduction_steps_file);

  for (UInt m = 0; m < nb_materials; ++m) {
    const Material & mat = model.getMaterial(m);
    if (mat.getName() == "gel")
      continue;
    // get the material parameters
    Real E = mat.getParam("E");
    Real max_damage = mat.getParam("max_damage");
    Real max_reductions = mat.getParam("max_reductions");

    // get the internal field that need to be set
    InternalField<UInt> & reduction_steps =
        const_cast<InternalField<UInt> &>(mat.getInternal<UInt>("damage_step"));
    InternalField<Real> & Sc =
        const_cast<InternalField<Real> &>(mat.getInternal<Real>("Sc"));

    InternalField<Real> & damage =
        const_cast<InternalField<Real> &>(mat.getInternal<Real>("damage"));
    // loop over all types in that material
    Real reduction_constant = mat.getParam("reduction_constant");
    for (auto element_type : reduction_steps.elementTypes(gt)) {
      const Array<UInt> & elem_filter = mat.getElementFilter(element_type, gt);
      if (!elem_filter.size())
        continue;
      Array<UInt> & reduction_step_array = reduction_steps(element_type, gt);
      Array<Real> & Sc_array = Sc(element_type, gt);
      Array<Real> & damage_array = damage(element_type, gt);
      const Array<Real> & D_tangent =
          mat.getInternal<Real>("tangent")(element_type, gt);
      const Array<Real> & eps_u =
          mat.getInternal<Real>("ultimate_strain")(element_type, gt);
      std::string line;
      for (UInt i = 0; i < reduction_step_array.size(); ++i) {

        std::getline(file_input, line);

        if (file_input.fail())
          AKANTU_EXCEPTION("Could not read data from file " +
                           reduction_steps_file);

        std::stringstream sstr(line);
        sstr >> reduction_step_array(i);
        if (reduction_step_array(i) == 0)
          continue;
        if (reduction_step_array(i) == max_reductions) {
          damage_array(i) = max_damage;
          Real previous_damage =
              1. -
              (1. / std::pow(reduction_constant, reduction_step_array(i) - 1));

          Sc_array(i) = eps_u(i) * (1. - previous_damage) * E * D_tangent(i) /
                        ((1. - previous_damage) * E + D_tangent(i));
        } else {
          damage_array(i) =
              1. - (1. / std::pow(reduction_constant, reduction_step_array(i)));
          Sc_array(i) = eps_u(i) * (1. - damage_array(i)) * E * D_tangent(i) /
                        ((1. - damage_array(i)) * E + D_tangent(i));
        }
      }
    }
  }
  // close the file
  file_input.close();
  return restart;
}
/* ------------------------------------------------------------------- */
Real RVETools::computePhaseVolume(const ID & mat_name) {
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  Real total_volume = 0;
  GhostType gt = _not_ghost;
  const Material & mat = model.getMaterial(mat_name);

  for (auto element_type : mesh.elementTypes(dim, gt, _ek_regular)) {
    const FEEngine & fe_engine = model.getFEEngine();
    const ElementTypeMapUInt & element_filter_map = mat.getElementFilter();
    if (!element_filter_map.exists(element_type, gt))
      continue;
    const Array<UInt> & elem_filter = mat.getElementFilter(element_type);
    if (!elem_filter.size())
      continue;
    Array<Real> volume(elem_filter.size() *
                           fe_engine.getNbIntegrationPoints(element_type),
                       1, 1.);

    total_volume += fe_engine.integrate(volume, element_type, gt, elem_filter);
  }

  // do not communicate if model is multi-scale
  if (not aka::is_of_type<SolidMechanicsModelRVE>(model)) {
    auto && comm = akantu::Communicator::getWorldCommunicator();
    comm.allReduce(total_volume, SynchronizerOperation::_sum);
  }

  return total_volume;
}
/* --------------------------------------------------------------------------
 */
void RVETools::applyEigenGradUinCracks(
    const Matrix<Real> prescribed_eigen_grad_u, const ID & material_name) {
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  GhostType gt = _not_ghost;
  Material & mat = model.getMaterial(material_name);
  const Real max_damage = mat.getParam("max_damage");
  const ElementTypeMapUInt & element_filter = mat.getElementFilter();
  for (auto element_type :
       element_filter.elementTypes(dim, gt, _ek_not_defined)) {

    if (!element_filter(element_type, gt).size())
      continue;
    const Array<Real> & damage =
        mat.getInternal<Real>("damage")(element_type, gt);
    const Real * dam_it = damage.storage();
    Array<Real> & eigen_gradu =
        mat.getInternal<Real>("eigen_grad_u")(element_type, gt);
    Array<Real>::matrix_iterator eigen_it = eigen_gradu.begin(dim, dim);
    Array<Real>::matrix_iterator eigen_end = eigen_gradu.end(dim, dim);
    for (; eigen_it != eigen_end; ++eigen_it, ++dam_it) {
      if (Math::are_float_equal(max_damage, *dam_it)) {
        Matrix<Real> & current_eigengradu = *eigen_it;
        current_eigengradu = prescribed_eigen_grad_u;
      }
    }
  }
}

/* --------------------------------------------------------------------------
 */
void RVETools::applyEigenGradUinCracks(
    const Matrix<Real> prescribed_eigen_grad_u,
    const ElementTypeMapUInt & critical_elements, bool to_add) {
  GhostType gt = _not_ghost;
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  for (auto element_type : mesh.elementTypes(dim, gt, _ek_not_defined)) {
    const Array<UInt> & critical_elements_vec =
        critical_elements(element_type, gt);
    const Array<UInt> & material_index_vec =
        model.getMaterialByElement(element_type, gt);
    const Array<UInt> & material_local_numbering_vec =
        model.getMaterialLocalNumbering(element_type, gt);

    for (UInt e = 0; e < critical_elements_vec.size(); ++e) {
      UInt element = critical_elements_vec(e);
      Material & material = model.getMaterial(material_index_vec(element));
      Array<Real> & eigen_gradu =
          material.getInternal<Real>("eigen_grad_u")(element_type, gt);
      UInt material_local_num = material_local_numbering_vec(element);
      Array<Real>::matrix_iterator eigen_it = eigen_gradu.begin(dim, dim);
      eigen_it += material_local_num;
      Matrix<Real> & current_eigengradu = *eigen_it;
      if (to_add)
        current_eigengradu += prescribed_eigen_grad_u;
      else
        current_eigengradu = prescribed_eigen_grad_u;
    }
  }
}

// code of Aurelia Cuba Ramos
// /*------------------------------------------------------------------------*/
//   void RVETools<dim>::fillCracks(ElementTypeMapReal & saved_damage) {
//   const UInt dim = model.getSpatialDimension();
//   const Material & mat_gel = model.getMaterial("gel");
//   const Real E_gel = mat_gel.getParam("E");
//   Real E_homogenized = 0.;
//   GhostType gt = _not_ghost;
//   const auto & mesh = model.getMesh();
//   for (UInt m = 0; m < model.getNbMaterials(); ++m) {
//     Material & mat = model.getMaterial(m);
//     if (mat.getName() == "gel")
//       continue;
//     const Real E = mat.getParam("E");
//     InternalField<Real> & damage = mat.getInternal<Real>("damage");
//     for (auto element_type : mesh.elementTypes(dim, gt, _ek_regular)) {
//       const Array<UInt> & elem_filter =
//       mat.getElementFilter(element_type, gt); if (!elem_filter.size())
//         continue;
//       Array<Real> & damage_vec = damage(element_type, gt);
//       Array<Real> & saved_damage_vec = saved_damage(element_type, gt);
//       for (UInt i = 0; i < damage_vec.size(); ++i) {
//         saved_damage_vec(elem_filter(i)) = damage_vec(i);
//         E_homogenized = (E_gel - E) * damage_vec(i) + E;
//         damage_vec(i) = 1. - (E_homogenized / E);
//       }
//     }
//   }
// }

// /*------------------------------------------------------------------------*/
// void RVETools<dim>::drainCracks(const ElementTypeMapReal &
// saved_damage) {
//   // model.dump();
//   const UInt dim = model.getSpatialDimension();
//   const auto & mesh = model.getMesh();
//   GhostType gt = _not_ghost;
//   for (UInt m = 0; m < model.getNbMaterials(); ++m) {
//     Material & mat = model.getMaterial(m);
//     if (mat.getName() == "gel")
//       continue;
//     else {
//       InternalField<Real> & damage = mat.getInternal<Real>("damage");
//       for (auto element_type : mesh.elementTypes(dim, gt,
//       _ek_not_defined)) {
//         const Array<UInt> & elem_filter =
//             mat.getElementFilter(element_type, gt);
//         if (!elem_filter.size())
//           continue;
//         Array<Real> & damage_vec = damage(element_type, gt);
//         const Array<Real> & saved_damage_vec =
//         saved_damage(element_type, gt); for (UInt i = 0; i <
//         damage_vec.size(); ++i) {
//           damage_vec(i) = saved_damage_vec(elem_filter(i));
//         }
//       }
//     }
//   }
//   // model.dump();
// }

/* --------------------------------------------------------------------------
 */
Real RVETools::computeSmallestElementSize() {

  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();

  // compute smallest element size
  const Array<Real> & pos = mesh.getNodes();
  Real el_h_min = std::numeric_limits<Real>::max();
  GhostType gt = _not_ghost;
  for (auto element_type : mesh.elementTypes(dim, gt)) {

    UInt nb_nodes_per_element = mesh.getNbNodesPerElement(element_type);
    UInt nb_element = mesh.getNbElement(element_type);

    Array<Real> X(0, nb_nodes_per_element * dim);
    model.getFEEngine().extractNodalToElementField(mesh, pos, X, element_type,
                                                   _not_ghost);

    Array<Real>::matrix_iterator X_el = X.begin(dim, nb_nodes_per_element);

    for (UInt el = 0; el < nb_element; ++el, ++X_el) {
      Real el_h = model.getFEEngine().getElementInradius(*X_el, element_type);
      if (el_h < el_h_min)
        el_h_min = el_h;
    }
  }

  if (not aka::is_of_type<SolidMechanicsModelRVE>(model)) {
    auto && comm = akantu::Communicator::getWorldCommunicator();
    comm.allReduce(el_h_min, SynchronizerOperation::_min);
  }

  return el_h_min;
}

/* ----------------------------------------------------------------------- */
void RVETools::computeChemExpansionLarive(
    const Real & delta_time_day, const Real & T, Real & expansion,
    const Real & eps_inf, const Real & time_ch_ref, const Real & time_lat_ref,
    const Real & U_C, const Real & U_L, const Real & T_ref) {
  AKANTU_DEBUG_IN();

  Real time_ch, time_lat, lambda, ksi, exp_ref;
  ksi = expansion / eps_inf;
  if (T == 0) {
    ksi += 0;
  } else {
    time_ch = time_ch_ref * std::exp(U_C * (1. / T - 1. / T_ref));
    time_lat = time_lat_ref * std::exp(U_L * (1. / T - 1. / T_ref));
    exp_ref = std::exp(-time_lat / time_ch);
    lambda = (1 + exp_ref) / (ksi + exp_ref);
    ksi += delta_time_day / time_ch * (1 - ksi) / lambda;
  }

  expansion = ksi * eps_inf;

  AKANTU_DEBUG_OUT();
}

/* ----------------------------------------------------------------------- */
void RVETools::computeChemExpansionArrhenius(const Real & delta_time_day,
                                             const Real & T, Real & expansion,
                                             const Real & k, const Real & Ea) {
  AKANTU_DEBUG_IN();

  if (T != 0) {
    const Real R = 8.3145; // J / mol / K
    expansion += delta_time_day * k * std::exp(-Ea / R / T);
  }

  AKANTU_DEBUG_OUT();
}
/* ----------------------------------------------------------------------- */
Real RVETools::computeDeltaEigenStrainThermal(const Real delta_time_day,
                                              const Real k,
                                              const Real activ_energy,
                                              const Real R, const Real T,
                                              Real & amount_reactive_particles,
                                              const Real saturation_const) {
  // compute increase in the eigen strain value for interval of time delta_time
  // as temperatures are stored in C, conversion to K is done

  Real delta_strain = amount_reactive_particles * k *
                      std::exp(-activ_energy / (R * T)) * delta_time_day;

  amount_reactive_particles -=
      std::exp(-activ_energy / (R * T)) * delta_time_day / saturation_const;

  if (amount_reactive_particles < 0.)
    amount_reactive_particles = 0.;
  return delta_strain;
}

/* -----------------------------------------------------------------------*/
Real RVETools::computeDeltaEigenStrainLinear(const Real delta_time,
                                             const Real k) {
  // compute increase in the eigen strain value for dt simply by deps = k *
  // delta_time

  Real delta_strain = k * delta_time;

  return delta_strain;
}

/* ---------------------------------------------------------------------- */
void RVETools::applyBoundaryConditionsRVE(
    const Matrix<Real> & displacement_gradient) {
  AKANTU_DEBUG_IN();

  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  // // transform into Green strain to exclude rotations
  // auto F = displacement_gradient + Matrix<Real>::eye(dim);
  // auto E = 0.5 * (F.transpose() * F - Matrix<Real>::eye(dim));
  // get the position of the nodes
  const Array<Real> & pos = mesh.getNodes();
  // storage for the coordinates of a given node and the displacement
  // that will be applied
  Vector<Real> x(dim);
  Vector<Real> appl_disp(dim);

  const auto & lower_bounds = mesh.getLowerBounds();
  for (auto node : this->corner_nodes) {
    x(0) = pos(node, 0);
    x(1) = pos(node, 1);
    x -= lower_bounds;
    appl_disp.mul<false>(displacement_gradient, x);
    (model.getBlockedDOFs())(node, 0) = true;
    (model.getDisplacement())(node, 0) = appl_disp(0);
    (model.getBlockedDOFs())(node, 1) = true;
    (model.getDisplacement())(node, 1) = appl_disp(1);
  }
  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
void RVETools::applyHomogeneousTemperature(const Real & temperature) {
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();

  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    Material & mat = model.getMaterial(m);

    for (auto el_type : mat.getElementFilter().elementTypes(dim)) {

      const auto & filter = mat.getElementFilter()(el_type);
      if (filter.size() == 0)
        continue;

      auto & deltas_T = mat.getArray<Real>("delta_T", el_type);

      for (auto && delta_T : deltas_T) {
        delta_T = temperature;
      }
    }
  }
}

/* --------------------------------------------------------------------------
 */
void RVETools::removeTemperature() {
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();

  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    Material & mat = model.getMaterial(m);

    for (auto el_type : mat.getElementFilter().elementTypes(dim)) {

      const auto & filter = mat.getElementFilter()(el_type);
      if (filter.size() == 0)
        continue;

      auto & deltas_T = mat.getArray<Real>("delta_T", el_type);

      for (auto && delta_T : deltas_T) {
        delta_T = 0.;
      }
    }
  }
}

/* --------------------------------------------------------------------------
 */
void RVETools::findCornerNodes() {
  AKANTU_DEBUG_IN();
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();

  // find corner nodes
  const auto & position = mesh.getNodes();
  const auto & lower_bounds = mesh.getLowerBounds();
  const auto & upper_bounds = mesh.getUpperBounds();

  switch (dim) {
  case 2: {
    corner_nodes.resize(4);
    corner_nodes.set(UInt(-1));
    for (auto && data : enumerate(make_view(position, dim))) {
      auto node = std::get<0>(data);
      const auto & X = std::get<1>(data);

      auto distance = X.distance(lower_bounds);
      // node 1 - left bottom corner
      if (Math::are_float_equal(distance, 0)) {
        corner_nodes(0) = node;
      }
      // node 2 - right bottom corner
      else if (Math::are_float_equal(X(_x), upper_bounds(_x)) &&
               Math::are_float_equal(X(_y), lower_bounds(_y))) {
        corner_nodes(1) = node;
      }
      // node 3 - right top
      else if (Math::are_float_equal(X(_x), upper_bounds(_x)) &&
               Math::are_float_equal(X(_y), upper_bounds(_y))) {
        corner_nodes(2) = node;
      }
      // node 4 - left top
      else if (Math::are_float_equal(X(_x), lower_bounds(_x)) &&
               Math::are_float_equal(X(_y), upper_bounds(_y))) {
        corner_nodes(3) = node;
      }
    }
    break;
  }
  case 3: {
    corner_nodes.resize(8);
    corner_nodes.set(UInt(-1));
    for (auto && data : enumerate(make_view(position, dim))) {
      auto node = std::get<0>(data);
      const auto & X = std::get<1>(data);

      auto distance = X.distance(lower_bounds);
      // node 1
      if (Math::are_float_equal(distance, 0)) {
        corner_nodes(0) = node;
      }
      // node 2
      else if (Math::are_float_equal(X(_x), upper_bounds(_x)) &&
               Math::are_float_equal(X(_y), lower_bounds(_y)) &&
               Math::are_float_equal(X(_z), lower_bounds(_z))) {
        corner_nodes(1) = node;
      }
      // node 3
      else if (Math::are_float_equal(X(_x), upper_bounds(_x)) &&
               Math::are_float_equal(X(_y), upper_bounds(_y)) &&
               Math::are_float_equal(X(_z), lower_bounds(_z))) {
        corner_nodes(2) = node;
      }
      // node 4
      else if (Math::are_float_equal(X(_x), lower_bounds(_x)) &&
               Math::are_float_equal(X(_y), upper_bounds(_y)) &&
               Math::are_float_equal(X(_z), lower_bounds(_z))) {
        corner_nodes(3) = node;
      }
      // node 5
      if (Math::are_float_equal(X(_x), lower_bounds(_x)) &&
          Math::are_float_equal(X(_y), lower_bounds(_y)) &&
          Math::are_float_equal(X(_z), upper_bounds(_z))) {
        corner_nodes(4) = node;
      }
      // node 6
      else if (Math::are_float_equal(X(_x), upper_bounds(_x)) &&
               Math::are_float_equal(X(_y), lower_bounds(_y)) &&
               Math::are_float_equal(X(_z), upper_bounds(_z))) {
        corner_nodes(5) = node;
      }
      // node 7
      else if (Math::are_float_equal(X(_x), upper_bounds(_x)) &&
               Math::are_float_equal(X(_y), upper_bounds(_y)) &&
               Math::are_float_equal(X(_z), upper_bounds(_z))) {
        corner_nodes(6) = node;
      }
      // node 8
      else if (Math::are_float_equal(X(_x), lower_bounds(_x)) &&
               Math::are_float_equal(X(_y), upper_bounds(_y)) &&
               Math::are_float_equal(X(_z), upper_bounds(_z))) {
        corner_nodes(7) = node;
      }
    }
    break;
  }
  }
  AKANTU_DEBUG_OUT();
}
/* --------------------------------------------------------------------------
 */
Real RVETools::averageScalarField(const ID & field_name) {
  AKANTU_DEBUG_IN();
  auto & fem = model.getFEEngine("SolidMechanicsFEEngine");
  Real average = 0;
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();

  GhostType gt = _not_ghost;
  for (auto element_type : mesh.elementTypes(dim, gt, _ek_not_defined)) {
    for (UInt m = 0; m < model.getNbMaterials(); ++m) {
      const auto & elem_filter =
          model.getMaterial(m).getElementFilter(element_type);
      if (!elem_filter.size())
        continue;
      const auto & scalar_field =
          model.getMaterial(m).getInternal<Real>(field_name)(element_type);
      Array<Real> int_scalar_vec(elem_filter.size(), 1, "int_of_scalar");

      fem.integrate(scalar_field, int_scalar_vec, 1, element_type, _not_ghost,
                    elem_filter);

      for (UInt k = 0; k < elem_filter.size(); ++k)
        average += int_scalar_vec(k);
    }
  }

  // compute total model volume
  if (!this->volume)
    computeModelVolume();

  return average / this->volume;

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
Real RVETools::averageTensorField(UInt row_index, UInt col_index,
                                  const ID & field_type) {
  AKANTU_DEBUG_IN();
  auto & fem = model.getFEEngine("SolidMechanicsFEEngine");
  Real average = 0;
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();

  GhostType gt = _not_ghost;
  for (auto element_type : mesh.elementTypes(dim, gt, _ek_not_defined)) {
    if (field_type == "stress") {
      for (UInt m = 0; m < model.getNbMaterials(); ++m) {
        const auto & stress_vec = model.getMaterial(m).getStress(element_type);
        const auto & elem_filter =
            model.getMaterial(m).getElementFilter(element_type);
        Array<Real> int_stress_vec(elem_filter.size(), dim * dim,
                                   "int_of_stress");

        fem.integrate(stress_vec, int_stress_vec, dim * dim, element_type,
                      _not_ghost, elem_filter);

        for (UInt k = 0; k < elem_filter.size(); ++k)
          average += int_stress_vec(k, row_index * dim + col_index);
      }
    } else if (field_type == "strain") {
      for (UInt m = 0; m < model.getNbMaterials(); ++m) {
        const auto & gradu_vec = model.getMaterial(m).getGradU(element_type);
        const auto & elem_filter =
            model.getMaterial(m).getElementFilter(element_type);
        Array<Real> int_gradu_vec(elem_filter.size(), dim * dim,
                                  "int_of_gradu");

        fem.integrate(gradu_vec, int_gradu_vec, dim * dim, element_type,
                      _not_ghost, elem_filter);

        for (UInt k = 0; k < elem_filter.size(); ++k)
          // averaging is done only for normal components, so stress and
          // strain are equal
          average += 0.5 * (int_gradu_vec(k, row_index * dim + col_index) +
                            int_gradu_vec(k, col_index * dim + row_index));
      }
    } else if (field_type == "eigen_grad_u") {
      for (UInt m = 0; m < model.getNbMaterials(); ++m) {
        const auto & eigen_gradu_vec = model.getMaterial(m).getInternal<Real>(
            "eigen_grad_u")(element_type);
        const auto & elem_filter =
            model.getMaterial(m).getElementFilter(element_type);
        Array<Real> int_eigen_gradu_vec(elem_filter.size(), dim * dim,
                                        "int_of_gradu");

        fem.integrate(eigen_gradu_vec, int_eigen_gradu_vec, dim * dim,
                      element_type, _not_ghost, elem_filter);

        for (UInt k = 0; k < elem_filter.size(); ++k)
          // averaging is done only for normal components, so stress and
          // strain are equal
          average += int_eigen_gradu_vec(k, row_index * dim + col_index);
      }
    } else {
      AKANTU_ERROR("Averaging not implemented for this field!!!");
    }
  }

  // compute total model volume
  if (!this->volume)
    computeModelVolume();

  return average / this->volume;

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
void RVETools::homogenizeStressField(Matrix<Real> & stress) {
  AKANTU_DEBUG_IN();
  stress(0, 0) = averageTensorField(0, 0, "stress");
  stress(1, 1) = averageTensorField(1, 1, "stress");
  stress(0, 1) = averageTensorField(0, 1, "stress");
  stress(1, 0) = averageTensorField(1, 0, "stress");
  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
void RVETools::setStiffHomogenDir(Matrix<Real> & stress) {
  AKANTU_DEBUG_IN();
  auto dim = stress.rows();
  Vector<Real> eigenvalues(dim);
  stress.eig(eigenvalues);
  Real hydrostatic_stress = 0;
  UInt denominator = 0;
  for (UInt i = 0; i < dim; ++i, ++denominator)
    hydrostatic_stress += eigenvalues(i);
  hydrostatic_stress /= denominator;
  AKANTU_DEBUG_OUT();
  this->tensile_homogenization = (hydrostatic_stress > 0);
}

/* --------------------------------------------------------------------------
 */
void RVETools::homogenizeStiffness(Matrix<Real> & C_macro, bool tensile_test) {
  AKANTU_DEBUG_IN();
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();

  AKANTU_DEBUG_ASSERT(dim == 2, "Is only implemented for 2D!!!");

  // apply three independent loading states to determine C
  // 1. eps_el = (0.001;0;0) 2. eps_el = (0,0.001,0) 3. eps_el = (0,0,0.001)

  // store and clear the eigenstrain
  // (to exclude stresses due to internal pressure)
  Array<Real> stored_eig(0, dim * dim, 0);
  storeEigenStrain(stored_eig);
  clearEigenStrain();

  // save nodal values before tests
  storeNodalFields();

  // storage for results of 3 different loading states
  UInt voigt_size = 1;
  switch (dim) {
  case 2: {
    voigt_size = VoigtHelper<2>::size;
    break;
  }
  case 3: {
    voigt_size = VoigtHelper<3>::size;
    break;
  }
  }
  Matrix<Real> stresses(voigt_size, voigt_size, 0.);
  Matrix<Real> strains(voigt_size, voigt_size, 0.);
  Matrix<Real> H(dim, dim, 0.);

  // indicate to material that its in the loading test
  UInt nb_materials = model.getNbMaterials();
  for (UInt m = 0; m < nb_materials; ++m) {
    Material & mat = model.getMaterial(m);
    if (aka::is_of_type<MaterialDamageIterativeOrthotropic<2>>(mat)) {
      mat.setParam("loading_test", true);
    }
  }

  // virtual test 1:
  H(0, 0) = 0.001 * (2 * tensile_test - 1);
  performVirtualTesting(H, stresses, strains, 0);

  // virtual test 2:
  H.zero();
  H(1, 1) = 0.001 * (2 * tensile_test - 1);
  performVirtualTesting(H, stresses, strains, 1);

  // virtual test 3:
  H.zero();
  H(0, 1) = 0.001;
  H(1, 0) = 0.001;
  performVirtualTesting(H, stresses, strains, 2);

  // indicate to material that its out of the loading test
  for (UInt m = 0; m < nb_materials; ++m) {
    Material & mat = model.getMaterial(m);
    if (aka::is_of_type<MaterialDamageIterativeOrthotropic<2>>(mat)) {
      mat.setParam("loading_test", false);
    }
  }

  // compute effective stiffness
  Matrix<Real> eps_inverse(voigt_size, voigt_size);
  eps_inverse.inverse(strains);

  // Make C matrix symmetric
  Matrix<Real> C_direct(voigt_size, voigt_size);
  C_direct.mul<false, false>(stresses, eps_inverse);
  for (UInt i = 0; i != voigt_size; ++i) {
    for (UInt j = 0; j != voigt_size; ++j) {
      C_macro(i, j) = 0.5 * (C_direct(i, j) + C_direct(j, i));
      C_macro(j, i) = C_macro(i, j);
    }
  }

  // return the nodal values and the eigenstrain
  restoreNodalFields();
  restoreEigenStrain(stored_eig);

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
void RVETools::performVirtualTesting(const Matrix<Real> & H,
                                     Matrix<Real> & eff_stresses,
                                     Matrix<Real> & eff_strains,
                                     const UInt test_no) {
  AKANTU_DEBUG_IN();
  auto & disp = model.getDisplacement();
  auto & boun = model.getBlockedDOFs();
  auto & ext_force = model.getExternalForce();
  disp.zero();
  boun.zero();
  ext_force.zero();

  applyBoundaryConditionsRVE(H);

  model.solveStep();

  // get average stress and strain
  eff_stresses(0, test_no) = averageTensorField(0, 0, "stress");
  eff_strains(0, test_no) = averageTensorField(0, 0, "strain");
  eff_stresses(1, test_no) = averageTensorField(1, 1, "stress");
  eff_strains(1, test_no) = averageTensorField(1, 1, "strain");
  eff_stresses(2, test_no) = averageTensorField(1, 0, "stress");
  eff_strains(2, test_no) = 2. * averageTensorField(1, 0, "strain");

  // restore historical internal fields
  restoreInternalFields();

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------- */
void RVETools::homogenizeEigenGradU(Matrix<Real> & eigen_gradu_macro) {
  AKANTU_DEBUG_IN();
  eigen_gradu_macro(0, 0) = averageTensorField(0, 0, "eigen_grad_u");
  eigen_gradu_macro(1, 1) = averageTensorField(1, 1, "eigen_grad_u");
  eigen_gradu_macro(0, 1) = averageTensorField(0, 1, "eigen_grad_u");
  eigen_gradu_macro(1, 0) = averageTensorField(1, 0, "eigen_grad_u");
  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
void RVETools::computeDamageRatio(Real & damage_ratio) {
  damage_ratio = 0.;
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    Material & mat = model.getMaterial(m);
    if (mat.getName() == "gel" || mat.getName() == "FE2_mat")
      continue;
    GhostType gt = _not_ghost;
    const ElementTypeMapArray<UInt> & filter_map = mat.getElementFilter();

    // Loop over the boundary element types
    for (auto & element_type : filter_map.elementTypes(dim, gt)) {
      const Array<UInt> & filter = filter_map(element_type);
      if (!filter_map.exists(element_type, gt))
        continue;
      if (filter.size() == 0)
        continue;

      const FEEngine & fe_engine = model.getFEEngine();
      auto & damage_array = mat.getInternal<Real>("damage")(element_type);

      damage_ratio +=
          fe_engine.integrate(damage_array, element_type, gt, filter);
    }
  }

  // do not communicate if model is multi-scale
  if (not aka::is_of_type<SolidMechanicsModelRVE>(model)) {
    auto && comm = akantu::Communicator::getWorldCommunicator();
    comm.allReduce(damage_ratio, SynchronizerOperation::_sum);
  }

  // compute total model volume
  if (!this->volume)
    computeModelVolume();

  damage_ratio /= this->volume;
}
/* --------------------------------------------------------------------------
 */
void RVETools::computeDamageRatioPerMaterial(Real & damage_ratio,
                                             const ID & material_name) {
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  GhostType gt = _not_ghost;
  Material & mat = model.getMaterial(material_name);
  const ElementTypeMapArray<UInt> & filter_map = mat.getElementFilter();
  const FEEngine & fe_engine = model.getFEEngine();
  damage_ratio = 0.;

  // Loop over the boundary element types
  for (auto & element_type : filter_map.elementTypes(dim, gt)) {
    const Array<UInt> & filter = filter_map(element_type);
    if (!filter_map.exists(element_type, gt))
      continue;
    if (filter.size() == 0)
      continue;

    auto & damage_array = mat.getInternal<Real>("damage")(element_type);

    damage_ratio += fe_engine.integrate(damage_array, element_type, gt, filter);
  }

  // do not communicate if model is multi-scale
  if (not aka::is_of_type<SolidMechanicsModelRVE>(model)) {
    auto && comm = akantu::Communicator::getWorldCommunicator();
    comm.allReduce(damage_ratio, SynchronizerOperation::_sum);
  }

  if (not this->phase_volumes.size())
    computePhaseVolumes();
  damage_ratio /= this->phase_volumes[material_name];
}

/* --------------------------------------------------------------------------
 */
void RVETools::computeCrackVolume(Real & crack_volume_ratio) {
  crack_volume_ratio = 0.;
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    Material & mat = model.getMaterial(m);
    if (mat.getName() == "gel" || mat.getName() == "FE2_mat")
      continue;
    GhostType gt = _not_ghost;
    const ElementTypeMapArray<UInt> & filter_map = mat.getElementFilter();

    // Loop over the boundary element types
    for (auto & element_type : filter_map.elementTypes(dim, gt)) {
      const Array<UInt> & filter = filter_map(element_type);
      if (!filter_map.exists(element_type, gt))
        continue;
      if (filter.size() == 0)
        continue;

      auto extra_volume_copy =
          mat.getInternal<Real>("extra_volume")(element_type);
      crack_volume_ratio += Math::reduce(extra_volume_copy);
    }
  }

  // do not communicate if model is multi-scale
  if (not aka::is_of_type<SolidMechanicsModelRVE>(model)) {
    auto && comm = akantu::Communicator::getWorldCommunicator();
    comm.allReduce(crack_volume_ratio, SynchronizerOperation::_sum);
  }

  // compute total model volume
  if (!this->volume)
    computeModelVolume();

  crack_volume_ratio /= this->volume;
}
/* --------------------------------------------------------------------------
 */
void RVETools::computeCrackVolumePerMaterial(Real & crack_volume,
                                             const ID & material_name) {
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  GhostType gt = _not_ghost;
  Material & mat = model.getMaterial(material_name);
  const ElementTypeMapArray<UInt> & filter_map = mat.getElementFilter();
  crack_volume = 0.;

  // Loop over the boundary element types
  for (auto & element_type : filter_map.elementTypes(dim, gt)) {
    const Array<UInt> & filter = filter_map(element_type);
    if (!filter_map.exists(element_type, gt))
      continue;
    if (filter.size() == 0)
      continue;

    auto extra_volume_copy =
        mat.getInternal<Real>("extra_volume")(element_type);

    crack_volume += Math::reduce(extra_volume_copy);
  }

  // do not communicate if model is multi-scale
  if (not aka::is_of_type<SolidMechanicsModelRVE>(model)) {
    auto && comm = akantu::Communicator::getWorldCommunicator();
    comm.allReduce(crack_volume, SynchronizerOperation::_sum);
  }

  if (not this->phase_volumes.size())
    computePhaseVolumes();
  crack_volume /= this->phase_volumes[material_name];
}

/* ------------------------------------------------------------------ */
void RVETools::dumpRVE() { model.dump(); }
/* ------------------------------------------------------------------ */
void RVETools::dumpRVE(UInt dump_nb) { model.dump(dump_nb); }

/* -------------------------------------------------------------------------
 */
void RVETools::communicateFlagsOnNodes() {
  auto & mesh = model.getMesh();
  auto & synch = mesh.getElementSynchronizer();
  NodesFlagUpdater nodes_flag_updater(mesh, synch,
                                      this->partition_border_nodes);
  nodes_flag_updater.fillPreventInsertion();
}
/* ------------------------------------------------------------------- */
void RVETools::communicateEffStressesOnNodes() {
  auto & mesh = model.getMesh();
  auto & synch = mesh.getElementSynchronizer();
  NodesEffStressUpdater nodes_eff_stress_updater(mesh, synch,
                                                 this->nodes_eff_stress);
  nodes_eff_stress_updater.updateMaxEffStressAtNodes();
}

/* ------------------------------------------------------------------ */
void RVETools::onNodesAdded(const Array<UInt> & new_nodes,
                            const NewNodesEvent &) {
  AKANTU_DEBUG_IN();
  if (new_nodes.size() == 0)
    return;
  // increase the internal arrays by the number of new_nodes
  UInt new_nb_nodes = this->modified_pos.size() + new_nodes.size();
  this->modified_pos.resize(new_nb_nodes);
  this->partition_border_nodes.resize(new_nb_nodes);
  this->nodes_eff_stress.resize(new_nb_nodes);
  this->crack_nodes.resize(new_nb_nodes);

  // function is activated only when expanding cohesive elements is on
  if (not this->cohesive_insertion)
    return;
  if (this->crack_central_nodes_ready)
    return;
  auto & mesh = model.getMesh();

  auto pos_it = make_view(mesh.getNodes(), mesh.getSpatialDimension()).begin();

  for (auto & new_node : new_nodes) {
    const Vector<Real> & new_node_coord = pos_it[new_node];
    auto central_nodes = this->crack_central_nodes;
    for (auto & central_node : central_nodes) {
      const Vector<Real> & central_node_coord = pos_it[central_node];
      if (new_node_coord == central_node_coord) {
        this->crack_central_nodes.push_back(new_node);
        break;
      }
    }
  }

  this->crack_central_nodes_ready = true;
  AKANTU_DEBUG_OUT();
}

/* ------------------------------------------------------------------------ */
void RVETools::updateElementGroup(const std::string group_name) {
  AKANTU_DEBUG_IN();
  auto & mesh = model.getMesh();
  AKANTU_DEBUG_ASSERT(mesh.elementGroupExists(group_name),
                      "Element group is not registered in the mesh");
  auto dim = mesh.getSpatialDimension();
  const GhostType gt = _not_ghost;
  auto && group = mesh.getElementGroup(group_name);
  auto && pos = mesh.getNodes();
  const auto pos_it = make_view(pos, dim).begin();

  for (auto & type : group.elementTypes(dim - 1)) {
    auto & facet_conn = mesh.getConnectivity(type, gt);
    const UInt nb_nodes_facet = facet_conn.getNbComponent();
    const auto facet_nodes_it = make_view(facet_conn, nb_nodes_facet).begin();
    AKANTU_DEBUG_ASSERT(
        type == _segment_2,
        "Currently update group works only for el type _segment_2");
    const auto element_ids = group.getElements(type, gt);
    for (auto && el_id : element_ids) {
      const auto connected_els = mesh.getElementToSubelement(type, gt)(el_id);
      for (auto && connected_el : connected_els) {

        auto type_solid = connected_el.type;
        auto & solid_conn = mesh.getConnectivity(type_solid, gt);
        const UInt nb_nodes_solid_el = solid_conn.getNbComponent();
        const auto solid_nodes_it =
            make_view(solid_conn, nb_nodes_solid_el).begin();
        Vector<UInt> facet_nodes = facet_nodes_it[el_id];
        Vector<UInt> solid_nodes = solid_nodes_it[connected_el.element];

        // check only the corner nodes of facets -central will not change
        for (auto f : arange(2)) {
          auto facet_node = facet_nodes(f);
          const Vector<Real> & facet_node_coords = pos_it[facet_node];
          for (auto s : arange(nb_nodes_solid_el)) {
            auto solid_node = solid_nodes(s);
            const Vector<Real> & solid_node_coords = pos_it[solid_node];
            if (solid_node_coords == facet_node_coords) {
              if (solid_node != facet_node) {
                // group.removeNode(facet_node);
                facet_conn(el_id, f) = solid_node;
                // group.addNode(solid_node, true);
              }
              break;
            }
          }
        }
      }
    }
  }

  group.optimize();
  AKANTU_DEBUG_OUT();
}

/* ------------------------------------------------------------------ */
void RVETools::applyEigenStrain(const Matrix<Real> & prestrain,
                                const ID & material_name) {
  AKANTU_DEBUG_IN();
  auto & mesh = model.getMesh();
  auto dim = mesh.getSpatialDimension();
  AKANTU_DEBUG_ASSERT(dim == 2, "This is 2D only!");

  // apply the new eigenstrain
  for (auto element_type :
       mesh.elementTypes(dim, _not_ghost, _ek_not_defined)) {
    Array<Real> & prestrain_vect = const_cast<Array<Real> &>(
        model.getMaterial(material_name)
            .getInternal<Real>("eigen_grad_u")(element_type));
    auto prestrain_it = prestrain_vect.begin(dim, dim);
    auto prestrain_end = prestrain_vect.end(dim, dim);

    for (; prestrain_it != prestrain_end; ++prestrain_it)
      (*prestrain_it) = prestrain;
  }
  AKANTU_DEBUG_OUT();
}

/* ------------------------------------------------------------------ */
void RVETools::clearEigenStrain(const ID & material_name) {
  AKANTU_DEBUG_IN();
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  Matrix<Real> zero_eigenstrain(dim, dim, 0.);
  GhostType gt = _not_ghost;
  for (auto element_type : mesh.elementTypes(dim, gt, _ek_not_defined)) {
    auto & prestrain_vect = const_cast<Array<Real> &>(
        model.getMaterial(material_name)
            .getInternal<Real>("eigen_grad_u")(element_type));
    auto prestrain_it = prestrain_vect.begin(dim, dim);
    auto prestrain_end = prestrain_vect.end(dim, dim);

    for (; prestrain_it != prestrain_end; ++prestrain_it)
      (*prestrain_it) = zero_eigenstrain;
  }
  AKANTU_DEBUG_OUT();
}
/* ------------------------------------------------------------------ */
void RVETools::storeEigenStrain(Array<Real> & stored_eig,
                                const ID & material_name) {
  AKANTU_DEBUG_IN();
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  GhostType gt = _not_ghost;
  UInt nb_quads = 0;
  for (auto element_type : mesh.elementTypes(dim, gt, _ek_not_defined)) {
    const UInt nb_gp = model.getFEEngine().getNbIntegrationPoints(element_type);
    nb_quads +=
        model.getMaterial(material_name).getElementFilter(element_type).size() *
        nb_gp;
  }
  stored_eig.resize(nb_quads);
  auto stored_eig_it = stored_eig.begin(dim, dim);
  for (auto element_type : mesh.elementTypes(dim, gt, _ek_not_defined)) {
    auto & prestrain_vect = const_cast<Array<Real> &>(
        model.getMaterial(material_name)
            .getInternal<Real>("eigen_grad_u")(element_type));
    auto prestrain_it = prestrain_vect.begin(dim, dim);
    auto prestrain_end = prestrain_vect.end(dim, dim);

    for (; prestrain_it != prestrain_end; ++prestrain_it, ++stored_eig_it)
      (*stored_eig_it) = (*prestrain_it);
  }
  AKANTU_DEBUG_OUT();
}
/* ------------------------------------------------------------------ */
void RVETools::restoreEigenStrain(Array<Real> & stored_eig,
                                  const ID & material_name) {
  AKANTU_DEBUG_IN();
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  GhostType gt = _not_ghost;
  auto stored_eig_it = stored_eig.begin(dim, dim);
  for (auto element_type : mesh.elementTypes(dim, gt, _ek_not_defined)) {
    auto & prestrain_vect = const_cast<Array<Real> &>(
        model.getMaterial(material_name)
            .getInternal<Real>("eigen_grad_u")(element_type));
    auto prestrain_it = prestrain_vect.begin(dim, dim);
    auto prestrain_end = prestrain_vect.end(dim, dim);

    for (; prestrain_it != prestrain_end; ++prestrain_it, ++stored_eig_it)
      (*prestrain_it) = (*stored_eig_it);
  }
  AKANTU_DEBUG_OUT();
}

/* ----------------------------------------------------------------- */
void RVETools::computeAveragePropertiesCohesiveModel(
    std::ofstream & file_output, Real time, Real offset) {
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();

  AKANTU_DEBUG_ASSERT(dim != 1, "Example does not work for 1D");

  Real av_strain_x = computeAverageStrain(_x, offset);
  Real av_strain_y = computeAverageStrain(_y, offset);

  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();

  if (dim == 2) {

    if (prank == 0)
      file_output << time << "," << av_strain_x << "," << av_strain_y << ","
                  << std::endl;
  }

  else {
    Real av_strain_z = computeAverageStrain(_z, offset);

    if (prank == 0)
      file_output << time << "," << av_strain_x << "," << av_strain_y << ","
                  << av_strain_z << std::endl;
  }
}
/* ------------------------------------------------------------------- */
void RVETools::insertCohesivesRandomly(const UInt & nb_insertions,
                                       std::string facet_mat_name,
                                       Real gap_ratio) {
  AKANTU_DEBUG_IN();

  // fill in the partition_border_nodes array
  communicateFlagsOnNodes();

  // pick central facets and neighbors
  pickFacetsRandomly(nb_insertions, facet_mat_name);

  insertOppositeFacetsAndCohesives();

  assignCrackNumbers();

  insertGap(gap_ratio);

  preventCohesiveInsertionInNeighbors();

  AKANTU_DEBUG_OUT();
}
/* ------------------------------------------------------------------- */
void RVETools::insertCohesivesRandomly3D(const UInt & nb_insertions,
                                         std::string facet_mat_name,
                                         Real /*gap_ratio*/) {
  AKANTU_DEBUG_IN();

  // fill in the partition_border_nodes array
  communicateFlagsOnNodes();

  // pick central facets and neighbors
  pickFacetsRandomly(nb_insertions, facet_mat_name);

  insertOppositeFacetsAndCohesives();

  assignCrackNumbers();

  // insertGap3D(gap_ratio);

  preventCohesiveInsertionInNeighbors();

  AKANTU_DEBUG_OUT();
}

/* ------------------------------------------------------------------- */
void RVETools::insertCohesiveLoops3D(const UInt & nb_insertions,
                                     std::string facet_mat_name,
                                     Real /*gap_ratio*/) {
  AKANTU_DEBUG_IN();

  // fill in the partition_border_nodes array
  communicateFlagsOnNodes();

  // pick central facets and neighbors
  auto inserted = closedFacetsLoopAroundPoint(nb_insertions, facet_mat_name);

  auto && comm = akantu::Communicator::getWorldCommunicator();
  comm.allReduce(inserted, SynchronizerOperation::_sum);
  AKANTU_DEBUG_ASSERT(inserted,
                      "No loading sites were inserted across the domain");

  insertOppositeFacetsAndCohesives();

  assignCrackNumbers();

  communicateCrackNumbers();

  // insertGap3D(gap_ratio);

  // preventCohesiveInsertionInNeighbors();

  AKANTU_DEBUG_OUT();
}

/* ------------------------------------------------------------------- */
template <UInt dim>
UInt RVETools::insertPreCracks(const UInt & nb_insertions,
                               std::string facet_mat_name) {
  AKANTU_DEBUG_IN();

  // fill in the partition_border_nodes array
  communicateFlagsOnNodes();

  // pick central facets and neighbors
  auto inserted = insertCohesiveLoops<dim>(nb_insertions, facet_mat_name);
  auto totally_inserted{inserted};
  auto && comm = akantu::Communicator::getWorldCommunicator();
  comm.allReduce(totally_inserted, SynchronizerOperation::_sum);
  AKANTU_DEBUG_ASSERT(totally_inserted,
                      "No loading sites were inserted across the domain");

  communicateCrackNumbers();

  AKANTU_DEBUG_OUT();
  return inserted;
}

template UInt RVETools::insertPreCracks<2>(const UInt & nb_insertions,
                                           std::string mat_name);
template UInt RVETools::insertPreCracks<3>(const UInt & nb_insertions,
                                           std::string mat_name);
/* -------------------------------------------------------------------------
 */
void RVETools::insertCohesivesByCoords(const Matrix<Real> & positions,
                                       Real gap_ratio) {
  AKANTU_DEBUG_IN();

  // fill in the border_nodes array
  communicateFlagsOnNodes();

  // pick the facets to duplicate
  pickFacetsByCoord(positions);

  insertOppositeFacetsAndCohesives();

  assignCrackNumbers();

  insertGap(gap_ratio);

  preventCohesiveInsertionInNeighbors();

  AKANTU_DEBUG_OUT();
}
/* ------------------------------------------------------------------- */
void RVETools::communicateCrackNumbers() {
  AKANTU_DEBUG_IN();

  if (model.getMesh().getCommunicator().getNbProc() == 1)
    return;

  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  const auto & coh_synch = coh_model.getCohesiveSynchronizer();
  CrackNumbersUpdater crack_numbers_updater(model, coh_synch);
  crack_numbers_updater.communicateCrackNumbers();
  AKANTU_DEBUG_OUT();
}
/* ----------------------------------------------------------------------- */
void RVETools::pickFacetsByCoord(const Matrix<Real> & positions) {
  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  auto & inserter = coh_model.getElementInserter();
  auto & mesh = model.getMesh();
  auto & mesh_facets = inserter.getMeshFacets();
  auto dim = mesh.getSpatialDimension();
  const GhostType gt = _not_ghost;
  auto type = *mesh.elementTypes(dim, gt, _ek_regular).begin();
  auto facet_type = Mesh::getFacetType(type);
  auto & doubled_facets =
      mesh_facets.createElementGroup("doubled_facets", dim - 1);
  auto & doubled_nodes = mesh.createNodeGroup("doubled_nodes");
  auto & facet_conn = mesh_facets.getConnectivity(facet_type, gt);
  auto & check_facets = inserter.getCheckFacets(facet_type, gt);
  const UInt nb_nodes_facet = facet_conn.getNbComponent();
  const auto facet_conn_it = make_view(facet_conn, nb_nodes_facet).begin();
  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();
  Vector<Real> bary_facet(dim);

  for (const auto & position : positions) {
    Real min_dist = std::numeric_limits<Real>::max();
    Element cent_facet;
    cent_facet.ghost_type = gt;
    for_each_element(
        mesh_facets,
        [&](auto && facet) {
          if (check_facets(facet.element)) {
            mesh_facets.getBarycenter(facet, bary_facet);
            auto dist = std::abs(bary_facet.distance(position));
            if (dist < min_dist) {
              min_dist = dist;
              cent_facet = facet;
            }
          }
        },
        _spatial_dimension = dim - 1, _ghost_type = _not_ghost);

    // communicate between processors for the element closest to the position
    auto local_min_dist = min_dist;
    auto && comm = akantu::Communicator::getWorldCommunicator();
    comm.allReduce(min_dist, SynchronizerOperation::_min);

    // let other processor insert the cohesive at this position
    if (local_min_dist != min_dist)
      continue;

    // eliminate possibility of inserting on a partition border
    bool border_facets{false};
    Vector<UInt> facet_nodes = facet_conn_it[cent_facet.element];
    for (auto node : facet_nodes) {
      if (this->partition_border_nodes(node))
        border_facets = true;
    }
    if (border_facets) {
      std::cout << "Facet at the position " << position(0) << "," << position(1)
                << " is close to the partition border. Skipping this position"
                << std::endl;
      continue;
    }

    // eliminate possibility of intersecting existing crack element
    bool intersection{false};
    for (auto node : facet_nodes) {
      if (this->crack_nodes(node))
        intersection = true;
    }
    if (intersection) {
      std::cout
          << "Facet at the position " << position(0) << "," << position(1)
          << " is intersecting another crack element. Skipping this position"
          << std::endl;
      continue;
    } else {
      auto success_with_neighbors = pickFacetNeighbors(cent_facet);

      if (success_with_neighbors) {
        doubled_facets.add(cent_facet);
        // add all facet nodes to the group
        for (auto node : arange(nb_nodes_facet)) {
          doubled_nodes.add(facet_conn(cent_facet.element, node));
          this->crack_nodes(facet_conn(cent_facet.element, node)) = true;
        }
        std::cout << "Proc " << prank << " placed 1 loading site" << std::endl;
        continue;
      } else
        continue;
    }
  }
}
/* ------------------------------------------------------------------- */
void RVETools::pickFacetsRandomly(UInt nb_insertions,
                                  std::string facet_mat_name) {
  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  auto & inserter = coh_model.getElementInserter();
  auto & mesh = model.getMesh();
  auto & mesh_facets = inserter.getMeshFacets();
  auto dim = mesh.getSpatialDimension();
  const GhostType gt = _not_ghost;
  auto type = *mesh.elementTypes(dim, gt, _ek_regular).begin();
  auto facet_type = Mesh::getFacetType(type);
  auto & doubled_facets =
      mesh_facets.createElementGroup("doubled_facets", dim - 1);
  auto & doubled_nodes = mesh.createNodeGroup("doubled_nodes");
  auto & matrix_elements =
      mesh_facets.createElementGroup("matrix_facets", dim - 1);
  auto facet_material_id = model.getMaterialIndex(facet_mat_name);
  auto & facet_conn = mesh_facets.getConnectivity(facet_type, gt);
  const UInt nb_nodes_facet = facet_conn.getNbComponent();
  const auto facet_conn_it = make_view(facet_conn, nb_nodes_facet).begin();
  auto & check_facets = inserter.getCheckFacets(facet_type, gt);
  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();

  if (not nb_insertions)
    return;

  Vector<Real> bary_facet(dim);
  for_each_element(
      mesh_facets,
      [&](auto && facet) {
        if (check_facets(facet.element)) {
          mesh_facets.getBarycenter(facet, bary_facet);
          auto & facet_material = coh_model.getFacetMaterial(
              facet.type, facet.ghost_type)(facet.element);
          if (facet_material == facet_material_id)
            matrix_elements.add(facet);
        }
      },
      _spatial_dimension = dim - 1, _ghost_type = _not_ghost);

  UInt nb_element = matrix_elements.getElements(facet_type).size();

  UInt seed = RandomGenerator<UInt>::seed();
  std::mt19937 random_generator(seed);
  std::uniform_int_distribution<> dis(0, nb_element - 1);

  if (not nb_element) {
    std::cout << "Proc " << prank << " couldn't place " << nb_insertions
              << " loading sites" << std::endl;
    return;
  }

  UInt already_inserted = 0;
  while (already_inserted < nb_insertions) {
    auto id = dis(random_generator);

    Element cent_facet;
    cent_facet.type = facet_type;
    cent_facet.element = matrix_elements.getElements(facet_type)(id);
    cent_facet.ghost_type = gt;

    // eliminate possibility of inserting on a partition border or
    // intersecting
    bool border_facets{false};
    Vector<UInt> facet_nodes = facet_conn_it[cent_facet.element];
    for (auto node : facet_nodes) {
      if (this->partition_border_nodes(node) or this->crack_nodes(node))
        border_facets = true;
    }
    if (border_facets)
      continue;
    else {
      auto success_with_neighbors = pickFacetNeighbors(cent_facet);

      if (success_with_neighbors) {
        already_inserted++;
        doubled_facets.add(cent_facet);
        // add all facet nodes to the group
        for (auto node : arange(nb_nodes_facet)) {
          doubled_nodes.add(facet_conn(cent_facet.element, node));
          this->crack_nodes(facet_conn(cent_facet.element, node)) = true;
        }
        std::cout << "Proc " << prank << " placed 1 loading site" << std::endl;
        continue;
      } else
        continue;
    }
  }
}
/* ------------------------------------------------------------------- */
UInt RVETools::closedFacetsLoopAroundPoint(UInt nb_insertions,
                                           std::string mat_name) {
  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  auto & inserter = coh_model.getElementInserter();
  auto & mesh = model.getMesh();
  auto & mesh_facets = inserter.getMeshFacets();
  auto dim = mesh.getSpatialDimension();
  const GhostType gt = _not_ghost;
  auto element_type = *mesh.elementTypes(dim, gt, _ek_regular).begin();
  auto facet_type = Mesh::getFacetType(element_type);
  auto & doubled_facets =
      mesh_facets.createElementGroup("doubled_facets", dim - 1);
  auto & doubled_nodes = mesh.createNodeGroup("doubled_nodes");
  auto material_id = model.getMaterialIndex(mat_name);
  auto & facet_conn = mesh_facets.getConnectivity(facet_type, gt);
  const UInt nb_nodes_facet = facet_conn.getNbComponent();
  const auto facet_conn_it = make_view(facet_conn, nb_nodes_facet).begin();
  auto & node_conn = mesh_facets.getConnectivity(_point_1, gt);
  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();

  const Real max_dot{0.7};

  if (not nb_insertions)
    return 0;

  CSR<Element> nodes_to_facets;
  MeshUtils::buildNode2Elements(mesh_facets, nodes_to_facets, dim - 1);

  std::set<UInt> matrix_nodes;

  for_each_element(
      mesh_facets,
      [&](auto && point_el) {
        // get the node number
        auto node = node_conn(point_el.element);
        // discard ghost nodes
        if (mesh.isPureGhostNode(node))
          goto nextnode;
        for (auto & facet : nodes_to_facets.getRow(node)) {
          auto & facet_material = coh_model.getFacetMaterial(
              facet.type, facet.ghost_type)(facet.element);
          if (facet_material != material_id) {
            goto nextnode;
          }
        }
        matrix_nodes.emplace(node);
      nextnode:;
      },
      _spatial_dimension = 0, _ghost_type = _not_ghost);

  UInt nb_nodes = matrix_nodes.size();
  UInt seed = RandomGenerator<UInt>::seed();
  std::mt19937 random_generator(seed);
  std::uniform_int_distribution<> dis(0, nb_nodes - 1);

  if (not nb_nodes) {
    std::cout << "Proc " << prank << " couldn't place " << nb_insertions
              << " loading sites" << std::endl;
    return 0;
  }

  UInt already_inserted = 0;
  std::set<UInt> left_nodes = matrix_nodes;
  while (already_inserted < nb_insertions) {

    if (not left_nodes.size()) {
      std::cout << "Proc " << prank << " inserted " << already_inserted
                << " loading sites out of " << nb_insertions << std::endl;
      return already_inserted;
    }
    auto id = dis(random_generator);
    auto it = matrix_nodes.begin();
    std::advance(it, id);

    UInt cent_node(*it);
    left_nodes.erase(*it);

    // not touching other cracks
    if (this->crack_nodes(cent_node))
      continue;

    // поехали
    CSR<Element> segments_to_nodes;
    MeshUtils::buildNode2Elements(mesh_facets, segments_to_nodes, dim - 2);
    auto && segments_to_node = segments_to_nodes.getRow(cent_node);
    Array<Element> segments_list;
    for (auto && segment : segments_to_node) {
      segments_list.push_back(segment);
    }

    // pick the first non-ghost segment in the list
    Element starting_segment;
    for (auto segment : segments_list) {
      if (segment.ghost_type == _not_ghost) {
        starting_segment = segment;
        break;
      }
    }

    // go to next node if starting segment could not be picked
    if (starting_segment == ElementNull) {
      std::cout
          << "Could not pick the starting segment, switching to the next node"
          << std::endl;
      continue;
    }

    auto & facets_to_segment =
        mesh_facets.getElementToSubelement(starting_segment);
    // pick the first good facet in the list
    Element starting_facet{ElementNull};
    for (const auto & facet : facets_to_segment) {
      // check if the facet and its nodes are ok
      if (isFacetAndNodesGood(facet, material_id, true)) {
        starting_facet = facet;
        break;
      }
    }

    // go to next node if starting facet could not be picked
    if (starting_facet == ElementNull) {
      std::cout
          << "Could not pick the starting facet, switching to the next node"
          << std::endl;
      continue;
    }

    // call the actual function that build the bridge
    auto facets_in_loop = findFacetsLoopFromSegment2Segment(
        starting_facet, starting_segment, starting_segment, cent_node, max_dot,
        material_id, true);

    // loop broke -> try another point
    if (not facets_in_loop.size()) {
      std::cout << "Couldn't close the facet loop, taking the next node"
                << std::endl;
      continue;
    }

    // otherwise champaigne
    for (auto & facet : facets_in_loop) {
      doubled_facets.add(facet);
      // add all facet nodes to the group
      for (auto node : arange(nb_nodes_facet)) {
        doubled_nodes.add(facet_conn(facet.element, node));
        this->crack_nodes(facet_conn(facet.element, node)) = true;
      }
    }
    this->crack_central_nodes.push_back(cent_node);

    already_inserted++;
  }
  std::cout << "Proc " << prank << " placed " << already_inserted
            << " loading site" << std::endl;
  return already_inserted;
}
/* ------------------------------------------------------------------- */
template <UInt dim>
UInt RVETools::insertCohesiveLoops(UInt nb_insertions, std::string mat_name) {
  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  auto & inserter = coh_model.getElementInserter();
  auto & mesh = model.getMesh();
  auto & mesh_facets = inserter.getMeshFacets();
  const GhostType gt = _not_ghost;
  auto element_type = *mesh.elementTypes(dim, gt, _ek_regular).begin();
  auto facet_type = Mesh::getFacetType(element_type);
  auto material_id = model.getMaterialIndex(mat_name);
  auto & facet_conn = mesh_facets.getConnectivity(facet_type, gt);
  const UInt nb_nodes_facet = facet_conn.getNbComponent();
  const auto facet_conn_it = make_view(facet_conn, nb_nodes_facet).begin();
  auto & node_conn = mesh_facets.getConnectivity(_point_1, gt);
  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();

  const Real max_dot{0.7};

  // instantiate crack_numbers array
  for (auto && type_facet : mesh_facets.elementTypes(dim - 1)) {
    ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);

    mesh.getDataPointer<UInt>("crack_numbers", type_cohesive);
  }
  ElementType type_cohesive = FEEngine::getCohesiveElementType(facet_type);
  auto & crack_numbers = mesh.getData<UInt>("crack_numbers", type_cohesive);

  CSR<Element> nodes_to_facets;
  MeshUtils::buildNode2Elements(mesh_facets, nodes_to_facets, dim - 1);

  std::set<UInt> matrix_nodes;

  for_each_element(
      mesh_facets,
      [&](auto && point_el) {
        // get the node number
        auto node = node_conn(point_el.element);
        // discard ghost nodes
        if (mesh.isPureGhostNode(node))
          goto nextnode;
        for (auto & facet : nodes_to_facets.getRow(node)) {
          auto & facet_material = coh_model.getFacetMaterial(
              facet.type, facet.ghost_type)(facet.element);
          if (facet_material != material_id) {
            goto nextnode;
          }
        }
        matrix_nodes.emplace(node);
      nextnode:;
      },
      _spatial_dimension = 0, _ghost_type = _not_ghost);

  UInt nb_nodes = matrix_nodes.size();
  UInt seed = RandomGenerator<UInt>::seed();
  std::mt19937 random_generator(seed);
  std::uniform_int_distribution<> dis(0, nb_nodes - 1);

  UInt already_inserted = 0;
  std::set<UInt> left_nodes = matrix_nodes;
  Array<Array<Element>> loops;
  while (already_inserted < nb_insertions and nb_nodes and left_nodes.size()) {

    auto id = dis(random_generator);
    auto it = matrix_nodes.begin();
    std::advance(it, id);

    UInt cent_node(*it);
    left_nodes.erase(*it);

    // not touching other cracks
    if (this->crack_nodes(cent_node))
      continue;

    // поехали
    CSR<Element> segments_to_nodes;
    MeshUtils::buildNode2Elements(mesh_facets, segments_to_nodes, dim - 2);
    auto && segments_to_node = segments_to_nodes.getRow(cent_node);
    Array<Element> segments_list;
    for (auto && segment : segments_to_node) {
      segments_list.push_back(segment);
    }

    // pick the first non-ghost segment in the list
    Element starting_segment;
    for (auto segment : segments_list) {
      if (segment.ghost_type == _not_ghost) {
        starting_segment = segment;
        break;
      }
    }

    // go to next node if starting segment could not be picked
    if (starting_segment == ElementNull) {
      std::cout
          << "Could not pick the starting segment, switching to the next node"
          << std::endl;
      continue;
    }

    auto & facets_to_segment =
        mesh_facets.getElementToSubelement(starting_segment);
    // pick the first good facet in the list
    Element starting_facet{ElementNull};
    for (const auto & facet : facets_to_segment) {
      // check if the facet and its nodes are ok
      if (isFacetAndNodesGood(facet, material_id, true)) {
        starting_facet = facet;
        break;
      }
    }

    // go to next node if starting facet could not be picked
    if (starting_facet == ElementNull) {
      std::cout
          << "Could not pick the starting facet, switching to the next node"
          << std::endl;
      continue;
    }

    // call the actual function that build the bridge
    auto facets_in_loop = findFacetsLoopFromSegment2Segment(
        starting_facet, starting_segment, starting_segment, cent_node, max_dot,
        material_id, true);

    // loop broke -> try another point
    if (not facets_in_loop.size()) {
      std::cout << "Couldn't close the facet loop, taking the next node"
                << std::endl;
      continue;
    }

    // otherwise champaigne
    loops.push_back(facets_in_loop);
    for (auto & facet : facets_in_loop) {
      // add all facet nodes to the group
      for (auto node : arange(nb_nodes_facet)) {
        this->crack_nodes(facet_conn(facet.element, node)) = true;
      }
    }
    this->crack_central_nodes.push_back(cent_node);

    already_inserted++;
  }

  // shift crack numbers by the number of cracks on lower processors
  auto num_cracks = loops.size();
  comm.exclusiveScan(num_cracks);

  // convert loops into a map
  std::map<UInt, UInt> facet_nbs_crack_nbs;
  for (auto && data : enumerate(loops)) {
    auto crack_nb = std::get<0>(data) + num_cracks;
    auto && facet_loop = std::get<1>(data);
    for (auto && facet : facet_loop) {
      facet_nbs_crack_nbs[facet.element] = crack_nb;
    }
  }

  // insert cohesives in a corresponding material
  auto & mat = model.getMaterial(material_id);
  auto * mat_coh = dynamic_cast<MaterialCohesiveLinearSequential<dim> *>(&mat);
  AKANTU_DEBUG_ASSERT(mat_coh != nullptr,
                      "Chosen material is not linear sequential");
  mat_coh->insertCohesiveElements(facet_nbs_crack_nbs, facet_type, false);
  // extend crack numbers array
  for (auto && pair : facet_nbs_crack_nbs) {
    auto && crack_nb = pair.second;
    crack_numbers.push_back(crack_nb);
  }

  inserter.insertElements();

  // fill up crack facets vector
  this->crack_facets.resize(loops.size());
  for (auto && data : zip(loops, this->crack_facets)) {
    auto && facet_loop = std::get<0>(data);
    auto && crack_facet_loop = std::get<1>(data);
    for (auto & facet : facet_loop) {
      Element cohesive{ElementNull};
      auto & connected_els = mesh_facets.getElementToSubelement(facet);
      for (auto & connected_el : connected_els) {
        if (mesh.getKind(connected_el.type) == _ek_cohesive) {
          cohesive = connected_el;
          break;
        }
      }
      AKANTU_DEBUG_ASSERT(cohesive != ElementNull,
                          "Connected cohesive element not identified");
      auto && connected_subels = mesh_facets.getSubelementToElement(cohesive);
      for (auto & connected_subel : connected_subels) {
        if (connected_subel.type != facet.type)
          continue;
        crack_facet_loop.push_back(connected_subel);
      }
    }
  }

  // split nodes into pairs and compute averaged normals
  identifyCrackCentralNodePairsAndNormals();

  std::cout << "Proc " << prank << " placed " << already_inserted
            << " loading site" << std::endl;
  return already_inserted;
}

template UInt RVETools::insertCohesiveLoops<2>(UInt nb_insertions,
                                               std::string mat_name);
template UInt RVETools::insertCohesiveLoops<3>(UInt nb_insertions,
                                               std::string mat_name);
/* ------------------------------------------------------------------- */
Array<Element> RVETools::findFacetsLoopFromSegment2Segment(
    Element directional_facet, Element starting_segment, Element ending_segment,
    UInt cent_node, Real max_dot, UInt material_id, bool check_crack_nodes) {

  auto & mesh = model.getMesh();
  auto & mesh_facets = mesh.getMeshFacets();
  auto dim = mesh.getSpatialDimension();
  AKANTU_DEBUG_ASSERT(
      dim == 3, "findFacetsLoopFromSegment2Segment currently supports only 3D");
  AKANTU_DEBUG_ASSERT(Mesh::getSpatialDimension(directional_facet.type) ==
                          dim - 1,
                      "Directional facet provided has wrong spatial dimension");

  if (Mesh::getSpatialDimension(starting_segment.type) != dim - 2 or
      Mesh::getSpatialDimension(ending_segment.type) != dim - 2) {
    AKANTU_EXCEPTION("Starting facet provided has wrong spatial dimension");
  }
  // get the point el corresponding to the central node
  CSR<Element> points_to_nodes;
  MeshUtils::buildNode2Elements(mesh_facets, points_to_nodes, 0);
  auto && points_to_node = points_to_nodes.getRow(cent_node);
  auto point_el = *points_to_node.begin();

  AKANTU_DEBUG_ASSERT(point_el.type == _point_1,
                      "Provided node element type is wrong");

  // check if the provided node is shared by two segments
  auto & segments_to_node = mesh_facets.getElementToSubelement(point_el);
  auto ret1 = std::find(segments_to_node.begin(), segments_to_node.end(),
                        starting_segment);
  AKANTU_DEBUG_ASSERT(ret1 != segments_to_node.end(),
                      "Starting segment doesn't contain the central node");
  auto ret2 = std::find(segments_to_node.begin(), segments_to_node.end(),
                        ending_segment);
  AKANTU_DEBUG_ASSERT(ret2 != segments_to_node.end(),
                      "Ending segment doesn't contain the central node");

  // check if the directional facet contain starting segment
  auto & facets_to_segment =
      mesh_facets.getElementToSubelement(starting_segment);
  auto ret3 = std::find(facets_to_segment.begin(), facets_to_segment.end(),
                        directional_facet);
  AKANTU_DEBUG_ASSERT(ret3 != facets_to_segment.end(),
                      "Starting facet doesn't contain the starting segment");

  // loop through facets to arrive from the starting segment to the endin
  bool loop_closed{false};
  auto current_facet(directional_facet);
  auto border_segment(starting_segment);
  Array<Element> facets_in_loop;
  Array<Element> segments_in_loop;
  while (not loop_closed and facets_in_loop.size() < 10) {

    // loop through all the neighboring facets and pick the prettiest
    Element next_facet(ElementNull);
    Real current_facet_indiam =
        MeshUtils::getInscribedCircleDiameter(model, current_facet);
    bool almost_closed{false};
    auto neighbor_facets = mesh_facets.getElementToSubelement(border_segment);
    for (auto neighbor_facet : neighbor_facets) {
      if (neighbor_facet == current_facet)
        continue;
      if (not isFacetAndNodesGood(neighbor_facet, material_id,
                                  check_crack_nodes))
        continue;
      // first check if it has the ending segment -> loop closed
      bool ending_segment_touched{false};
      if (facets_in_loop.size()) {
        const Vector<Element> neighbors_segments =
            mesh_facets.getSubelementToElement(neighbor_facet);
        auto ret = std::find(neighbors_segments.begin(),
                             neighbors_segments.end(), ending_segment);
        if (ret != neighbors_segments.end()) {
          ending_segment_touched = true;
          almost_closed = true;
        }
      }

      // conditions on the angle between facets
      Real neighbor_facet_indiam =
          MeshUtils::getInscribedCircleDiameter(model, neighbor_facet);
      Real hypotenuse = sqrt(neighbor_facet_indiam * neighbor_facet_indiam +
                             current_facet_indiam * current_facet_indiam) /
                        2;
      auto dist = MeshUtils::distanceBetweenTwoFacetsAlligned(
          mesh_facets, current_facet, neighbor_facet);
      if (dist < hypotenuse) {
        continue;
      }

      // get abs of dot product between two normals
      Real dot = MeshUtils::cosSharpAngleBetween2Facets(model, current_facet,
                                                        neighbor_facet);
      if (dot > max_dot) {
        next_facet = neighbor_facet;
        if (ending_segment_touched) {
          facets_in_loop.push_back(neighbor_facet);
          loop_closed = true;
          break;
        }
      }
    }

    if (next_facet == ElementNull) {
      facets_in_loop.clear();
      return facets_in_loop;
    }

    if (almost_closed & not loop_closed) {
      facets_in_loop.clear();
      return facets_in_loop;
    }

    if (loop_closed) {
      return facets_in_loop;
    }

    // otherwise make this facet a current one
    current_facet = next_facet;
    facets_in_loop.push_back(current_facet);
    segments_in_loop.push_back(border_segment);

    // identify the next border segment
    auto previous_border_segment = border_segment;
    border_segment = ElementNull;
    const Vector<Element> current_segments =
        mesh_facets.getSubelementToElement(current_facet);
    for (auto & segment : current_segments) {
      if (segment == previous_border_segment)
        continue;
      // check if the segment includes the central node
      auto && nodes_to_segment = mesh_facets.getConnectivity(segment);
      bool includes_cent_node{false};
      for (auto & node : nodes_to_segment) {
        if (node == cent_node) {
          includes_cent_node = true;
          break;
        }
      }
      if (not includes_cent_node)
        continue;

      // check if the segment was previously added
      if (segments_in_loop.find(segment) != UInt(-1))
        continue;

      border_segment = segment;
      break;
    }
    if (border_segment == ElementNull) {
      facets_in_loop.clear();
      return facets_in_loop;
    }
  }
  return facets_in_loop;
}
/* ------------------------------------------------------------------- */
Array<Element>
RVETools::findFacetsLoopByGraphByDist(const Array<Element> & limit_facets,
                                      const Array<Element> & limit_segments,
                                      const UInt & cent_node) {

  auto & mesh = model.getMesh();
  auto & mesh_facets = mesh.getMeshFacets();
  auto dim = mesh.getSpatialDimension();
  Array<Element> facets_in_loop;
  auto && starting_facet = limit_facets(0);
  auto && ending_facet = limit_facets(1);
  auto && starting_segment = limit_segments(0);
  auto && ending_segment = limit_segments(1);
  Real max_dist = std::numeric_limits<Real>::min();

  // Create a struct to hold properties for each vertex
  typedef struct vertex_properties {
    Element segment;
  } vertex_properties_t;

  // Create a struct to hold properties for each edge
  typedef struct edge_properties {
    Element facet;
    Real weight;
  } edge_properties_t;

  // Define the type of the graph
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                                vertex_properties_t, edge_properties_t>
      graph_t;
  typedef graph_t::vertex_descriptor vertex_descriptor_t;
  typedef graph_t::edge_descriptor edge_descriptor_t;
  typedef boost::graph_traits<graph_t>::edge_iterator edge_iter;
  typedef boost::graph_traits<graph_t>::vertex_iterator vertex_iter;
  graph_t g;

  // prepare the list of facets sharing the central node
  CSR<Element> facets_to_nodes;
  MeshUtils::buildNode2Elements(mesh_facets, facets_to_nodes, dim - 1);
  auto && facets_to_node = facets_to_nodes.getRow(cent_node);
  std::set<Element> facets_added;
  std::map<Element, vertex_descriptor_t> segment_vertex_map;

  // get the point el corresponding to the central node
  CSR<Element> points_to_nodes;
  MeshUtils::buildNode2Elements(mesh_facets, points_to_nodes, 0);
  auto && points_to_node = points_to_nodes.getRow(cent_node);
  auto point_el = *points_to_node.begin();
  // get list of all segments connected to the node
  auto && segm_to_point = mesh_facets.getElementToSubelement(point_el);

  // add starting and ending segment to the graph
  vertex_descriptor_t v_start = boost::add_vertex(g);
  g[v_start].segment = starting_segment;
  vertex_descriptor_t v_end = boost::add_vertex(g);
  g[v_end].segment = ending_segment;
  segment_vertex_map[starting_segment] = v_start;
  segment_vertex_map[ending_segment] = v_end;

  for (auto && facet : facets_to_node) {
    // discard undesired facets first
    if (facet == starting_facet or facet == ending_facet)
      continue;
    if (not isFacetAndNodesGood(facet))
      continue;
    if (belong2SameElement(starting_facet, facet) or
        belong2SameElement(ending_facet, facet))
      continue;

    // identify 2 side facets
    Array<Element> side_segments;
    auto && facet_segments = mesh_facets.getSubelementToElement(facet);
    for (auto & segment : facet_segments) {
      auto ret = std::find(segm_to_point.begin(), segm_to_point.end(), segment);
      if (ret != segm_to_point.end()) {
        side_segments.push_back(segment);
      }
    }
    AKANTU_DEBUG_ASSERT(side_segments.size() == 2,
                        "Number of side facets is wrong");

    // add corresponding vertex (if not there)and edge into graph
    Array<vertex_descriptor_t> v(2);
    for (UInt i = 0; i < 2; i++) {
      auto it = segment_vertex_map.find(side_segments(i));
      if (it != segment_vertex_map.end()) {
        v(i) = it->second;
      } else {
        v(i) = boost::add_vertex(g);
        g[v(i)].segment = side_segments(i);
        segment_vertex_map[side_segments(i)] = v(i);
      }
    }

    // add corresponding edge into the graph
    auto ret = facets_added.emplace(facet);
    if (ret.second) {
      // compute facet's area to use as a weight
      auto dist = MeshUtils::distanceBetweenTwoFacets(mesh_facets,
                                                      starting_facet, facet);
      dist +=
          MeshUtils::distanceBetweenTwoFacets(mesh_facets, ending_facet, facet);
      if (dist > max_dist)
        max_dist = dist;

      // add edge
      std::pair<edge_descriptor_t, bool> e = boost::add_edge(v(0), v(1), g);
      g[e.first].weight = dist;
      // g[e.first].weight = facet_area;
      g[e.first].facet = facet;
    }
  }

  // normalize and flip weights
  edge_iter ei, ei_end;
  for (std::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
    auto weight = g[*ei].weight;
    weight = max_dist - weight;
    auto & fe_engine = model.getFEEngine("FacetsFEEngine");
    auto facet_area = MeshUtils::getFacetArea(fe_engine, g[*ei].facet);
    g[*ei].weight = weight + sqrt(facet_area);
  }

  // Print out some useful information
  std::cout << "Starting segm " << starting_segment.element << " ending segm "
            << ending_segment.element << std::endl;
  std::cout << "Graph:" << std::endl;
  for (std::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
    std::cout << g[*ei].facet.element << " ";
  }
  std::cout << std::endl;
  std::cout << "num_verts: " << boost::num_vertices(g) << std::endl;
  std::cout << "num_edges: " << boost::num_edges(g) << std::endl;

  // BGL Dijkstra's Shortest Paths here...
  std::vector<Real> distances(boost::num_vertices(g));
  std::vector<vertex_descriptor_t> predecessors(boost::num_vertices(g));

  boost::dijkstra_shortest_paths(
      g, v_start,
      boost::weight_map(boost::get(&edge_properties_t::weight, g))
          .distance_map(boost::make_iterator_property_map(
              distances.begin(), boost::get(boost::vertex_index, g)))
          .predecessor_map(boost::make_iterator_property_map(
              predecessors.begin(), boost::get(boost::vertex_index, g))));
  std::cout << "distances and parents:" << std::endl;
  vertex_iter vi, vend;
  for (boost::tie(vi, vend) = boost::vertices(g); vi != vend; ++vi) {
    std::cout << "distance(" << g[*vi].segment.element
              << ") = " << distances[*vi] << ", ";
    std::cout << "parent = " << g[predecessors[*vi]].segment.element
              << std::endl;
  }

  // Extract the shortest path from v1 to v3.
  typedef std::vector<edge_descriptor_t> path_t;
  path_t path;

  vertex_descriptor_t v = v_end;
  for (vertex_descriptor_t u = predecessors[v]; u != v;
       v = u, u = predecessors[v]) {
    std::pair<edge_descriptor_t, bool> edge_pair = boost::edge(u, v, g);
    path.push_back(edge_pair.first);
  }

  std::cout << std::endl;
  std::cout << "Shortest Path from segment " << starting_segment.element
            << " to " << ending_segment.element << ":" << std::endl;
  for (path_t::reverse_iterator riter = path.rbegin(); riter != path.rend();
       ++riter) {
    vertex_descriptor_t u_tmp = boost::source(*riter, g);
    vertex_descriptor_t v_tmp = boost::target(*riter, g);
    edge_descriptor_t e_tmp = boost::edge(u_tmp, v_tmp, g).first;

    std::cout << "  " << g[u_tmp].segment.element << " -> "
              << g[v_tmp].segment.element << " through "
              << g[e_tmp].facet.element << "    (area: " << g[e_tmp].weight
              << ")" << std::endl;
    facets_in_loop.push_back(g[e_tmp].facet);
  }
  return facets_in_loop;
}
/* ------------------------------------------------------------------- */
Array<Element> RVETools::findFacetsLoopByGraphByArea(
    const Array<Element> & limit_facets, const Array<Element> & limit_segments,
    const Array<Element> & preinserted_facets, const UInt & cent_node) {

  auto & mesh = model.getMesh();
  auto & mesh_facets = mesh.getMeshFacets();
  auto dim = mesh.getSpatialDimension();
  Array<Element> facets_in_loop;

  // find two best neighbors
  Array<Element> best_neighbors(2, 1, ElementNull);
  Array<Element> neighbors_borders(2, 1, ElementNull);
  for (UInt i = 0; i < 2; i++) {
    Element best_neighbor{ElementNull};
    // if the facet is already provided - use it
    if (preinserted_facets(i) != ElementNull) {
      best_neighbor = preinserted_facets(i);
    } else {
      // if not - search for it
      auto limit_facet_indiam =
          MeshUtils::getInscribedCircleDiameter(model, limit_facets(i));
      auto max_dot = std::numeric_limits<Real>::min();
      auto && neighbor_facets =
          mesh_facets.getElementToSubelement(limit_segments(i));
      for (auto && neighbor_facet : neighbor_facets) {
        if (neighbor_facet == limit_facets(0) or
            neighbor_facet == limit_facets(1))
          continue;
        if (not isFacetAndNodesGood(neighbor_facet))
          continue;
        // conditions on the angle between facets
        Real neighbor_facet_indiam =
            MeshUtils::getInscribedCircleDiameter(model, neighbor_facet);
        Real hypotenuse = sqrt(neighbor_facet_indiam * neighbor_facet_indiam +
                               limit_facet_indiam * limit_facet_indiam) /
                          2;
        auto dist = MeshUtils::distanceBetweenTwoFacetsAlligned(
            mesh_facets, limit_facets(i), neighbor_facet);
        if (dist < hypotenuse)
          continue;

        // find neighbor with the biggest thing
        // get abs of dot product between two normals
        Real dot = MeshUtils::cosSharpAngleBetween2Facets(
            model, limit_facets(i), neighbor_facet);
        if (dot > max_dot) {
          max_dot = dot;
          best_neighbor = neighbor_facet;
        }
      }
    }
    if (best_neighbor == ElementNull) {
      return facets_in_loop;
    } else {
      best_neighbors(i) = best_neighbor;
    }

    // find new limiting segments
    auto && neighbor_segments =
        mesh_facets.getSubelementToElement(best_neighbor);
    for (auto & segment : neighbor_segments) {
      if (segment == limit_segments(i))
        continue;

      // check if the segment includes the central node
      auto && nodes_to_segment = mesh_facets.getConnectivity(segment);
      bool includes_cent_node{false};
      for (auto & node : nodes_to_segment) {
        if (node == cent_node) {
          includes_cent_node = true;
          break;
        }
      }
      if (not includes_cent_node)
        continue;

      neighbors_borders(i) = segment;
      break;
    }
  }
  // check if borders are not the same segment
  if (neighbors_borders(0) == neighbors_borders(1)) {
    for (auto && best_neighbor : best_neighbors) {
      facets_in_loop.push_back(best_neighbor);
    }
  }

  // Create a struct to hold properties for each vertex
  typedef struct vertex_properties {
    Element segment;
  } vertex_properties_t;

  // Create a struct to hold properties for each edge
  typedef struct edge_properties {
    Element facet;
    Real weight;
  } edge_properties_t;

  // Define the type of the graph
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                                vertex_properties_t, edge_properties_t>
      graph_t;
  typedef graph_t::vertex_descriptor vertex_descriptor_t;
  typedef graph_t::edge_descriptor edge_descriptor_t;
  graph_t g;

  // prepare the list of facets sharing the central node
  CSR<Element> facets_to_nodes;
  MeshUtils::buildNode2Elements(mesh_facets, facets_to_nodes, dim - 1);
  auto && facets_to_node = facets_to_nodes.getRow(cent_node);
  std::set<Element> facets_added;
  std::map<Element, vertex_descriptor_t> segment_vertex_map;

  // get the point el corresponding to the central node
  CSR<Element> points_to_nodes;
  MeshUtils::buildNode2Elements(mesh_facets, points_to_nodes, 0);
  auto && points_to_node = points_to_nodes.getRow(cent_node);
  auto point_el = *points_to_node.begin();
  // get list of all segments connected to the node
  auto && segm_to_point = mesh_facets.getElementToSubelement(point_el);

  // add starting and ending segment to the graph
  vertex_descriptor_t v_start = boost::add_vertex(g);
  g[v_start].segment = neighbors_borders(0);
  vertex_descriptor_t v_end = boost::add_vertex(g);
  g[v_end].segment = neighbors_borders(1);
  segment_vertex_map[neighbors_borders(0)] = v_start;
  segment_vertex_map[neighbors_borders(1)] = v_end;

  // compute inscribed circles diameters
  Array<Real> best_neighbors_indiam;
  for (auto && best_neighbor : best_neighbors) {
    auto neighbor_facet_indiam =
        MeshUtils::getInscribedCircleDiameter(model, best_neighbor);
    best_neighbors_indiam.push_back(neighbor_facet_indiam);
  }

  for (auto && facet : facets_to_node) {
    // discard undesired facets first
    if (facet == best_neighbors(0) or facet == best_neighbors(1) or
        facet == limit_facets(0) or facet == limit_facets(1))
      continue;
    if (not isFacetAndNodesGood(facet))
      continue;
    if (belong2SameElement(best_neighbors(0), facet) or
        belong2SameElement(best_neighbors(1), facet) or
        belong2SameElement(limit_facets(0), facet) or
        belong2SameElement(limit_facets(1), facet))
      continue;

    // discard short distances to starting and ending facets
    bool short_dist{false};
    auto facet_indiam = MeshUtils::getInscribedCircleDiameter(model, facet);
    for (UInt i = 0; i < 2; i++) {
      auto hypotenuse =
          sqrt(best_neighbors_indiam(i) * best_neighbors_indiam(i) +
               facet_indiam * facet_indiam) /
          2;
      auto dist = MeshUtils::distanceBetweenTwoFacets(mesh_facets,
                                                      best_neighbors(i), facet);

      if (dist <= hypotenuse) {
        short_dist = true;
        break;
      }
    }
    if (short_dist)
      continue;

    // identify 2 side facets
    Array<Element> side_segments;
    auto && facet_segments = mesh_facets.getSubelementToElement(facet);
    for (auto & segment : facet_segments) {
      auto ret = std::find(segm_to_point.begin(), segm_to_point.end(), segment);
      if (ret != segm_to_point.end()) {
        side_segments.push_back(segment);
      }
    }
    AKANTU_DEBUG_ASSERT(side_segments.size() == 2,
                        "Number of side facets is wrong");

    // add corresponding vertex (if not there)and edge into graph
    Array<vertex_descriptor_t> v(2);
    for (UInt i = 0; i < 2; i++) {
      auto it = segment_vertex_map.find(side_segments(i));
      if (it != segment_vertex_map.end()) {
        v(i) = it->second;
      } else {
        v(i) = boost::add_vertex(g);
        g[v(i)].segment = side_segments(i);
        segment_vertex_map[side_segments(i)] = v(i);
      }
    }

    // add corresponding edge into the graph
    auto ret = facets_added.emplace(facet);
    if (ret.second) {
      // compute facet's area to use as a weight
      auto & fe_engine = model.getFEEngine("FacetsFEEngine");
      auto facet_area = MeshUtils::getFacetArea(fe_engine, facet);

      // add edge
      std::pair<edge_descriptor_t, bool> e = boost::add_edge(v(0), v(1), g);
      g[e.first].weight = facet_area;
      g[e.first].facet = facet;
    }
  }

  // // Print out some useful information
  // std::cout << "Starting facet " << best_neighbors(0).element
  //           << " ending facet " << best_neighbors(0).element << std::endl;
  // std::cout << "Graph:" << std::endl;
  // // boost::print_graph(g, boost::get(&vertex_properties_t::segment, g));
  // edge_iter ei, ei_end;
  // for (std::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
  //   std::cout << g[*ei].facet.element << ",";
  // }
  // std::cout << std::endl;
  // for (std::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
  //   vertex_descriptor_t u = boost::source(*ei, g);
  //   vertex_descriptor_t v = boost::target(*ei, g);
  //   std::cout << g[*ei].facet.element << " from " << g[u].segment.element
  //             << "->" << g[v].segment.element << " area " << g[*ei].weight
  //             << std::endl;
  // }
  // std::cout << "num_verts: " << boost::num_vertices(g) << std::endl;
  // std::cout << "num_edges: " << boost::num_edges(g) << std::endl;
  // std::cout << std::endl;

  // BGL Dijkstra's Shortest Paths here...
  std::vector<Real> distances(boost::num_vertices(g));
  std::vector<vertex_descriptor_t> predecessors(boost::num_vertices(g));

  boost::dijkstra_shortest_paths(
      g, v_end,
      boost::weight_map(boost::get(&edge_properties_t::weight, g))
          .distance_map(boost::make_iterator_property_map(
              distances.begin(), boost::get(boost::vertex_index, g)))
          .predecessor_map(boost::make_iterator_property_map(
              predecessors.begin(), boost::get(boost::vertex_index, g))));
  // Extract the shortest path from v1 to v3.
  typedef std::vector<edge_descriptor_t> path_t;
  path_t path;

  vertex_descriptor_t v = v_start;
  for (vertex_descriptor_t u = predecessors[v]; u != v;
       v = u, u = predecessors[v]) {
    std::pair<edge_descriptor_t, bool> edge_pair = boost::edge(u, v, g);
    path.push_back(edge_pair.first);
  }

  if (path.size()) {
    // add 2 determenistic facets
    for (auto && best_neighbor : best_neighbors) {
      facets_in_loop.push_back(best_neighbor);
    }
    // add the shortest path facets
    for (path_t::reverse_iterator riter = path.rbegin(); riter != path.rend();
         ++riter) {
      vertex_descriptor_t u_tmp = boost::source(*riter, g);
      vertex_descriptor_t v_tmp = boost::target(*riter, g);
      edge_descriptor_t e_tmp = boost::edge(u_tmp, v_tmp, g).first;

      facets_in_loop.push_back(g[e_tmp].facet);
    }
  }
  return facets_in_loop;
}
/* ------------------------------------------------------------------- */
Array<Element> RVETools::findStressedFacetLoopAroundNode(
    const Array<Element> & limit_facets, const Array<Element> & limit_segments,
    const UInt & cent_node, Real min_dot) {

  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  auto & mesh = model.getMesh();
  auto & mesh_facets = mesh.getMeshFacets();
  auto dim = mesh.getSpatialDimension();
  AKANTU_DEBUG_ASSERT(limit_facets(0) != limit_facets(1),
                      "Limit facets are equal");
  AKANTU_DEBUG_ASSERT(limit_segments(0) != limit_segments(1),
                      "Limit segments are equal");
  AKANTU_DEBUG_ASSERT(
      dim == 3, "findFacetsLoopFromSegment2Segment currently supports only 3D");
  // get the point el corresponding to the central node
  CSR<Element> points_to_nodes;
  MeshUtils::buildNode2Elements(mesh_facets, points_to_nodes, 0);
  auto && points_to_node = points_to_nodes.getRow(cent_node);
  auto point_el = *points_to_node.begin();
  auto & segments_to_point = mesh_facets.getElementToSubelement(point_el);

  AKANTU_DEBUG_ASSERT(point_el.type == _point_1,
                      "Provided node element type is wrong");

  for (UInt i = 0; i < 2; i++) {
    // check the dimensions
    AKANTU_DEBUG_ASSERT(
        Mesh::getSpatialDimension(limit_facets(i).type) == dim - 1,
        "Facet" << limit_facets(i) << " has wrong spatial dimension");
    AKANTU_DEBUG_ASSERT(
        Mesh::getSpatialDimension(limit_segments(i).type) == dim - 2,
        "Segment" << limit_segments(i) << " has wrong spatial dimension");

    // check if the provided node is shared by two segments
    auto ret = std::find(segments_to_point.begin(), segments_to_point.end(),
                         limit_segments(i));
    AKANTU_DEBUG_ASSERT(ret != segments_to_point.end(),
                        "Segment " << limit_segments(i)
                                   << " doesn't contain the central node");

    // check if limit facets contain limit segments
    auto & facets_to_segment =
        mesh_facets.getElementToSubelement(limit_segments(i));
    auto ret1 = std::find(facets_to_segment.begin(), facets_to_segment.end(),
                          limit_facets(i));
    AKANTU_DEBUG_ASSERT(ret1 != facets_to_segment.end(),
                        "Facet " << limit_facets(i)
                                 << " doesn't contain the segment "
                                 << limit_segments(i));
  }
  Array<Element> facets_in_loop;

  // Create a struct to hold properties for each vertex
  typedef struct vertex_properties {
    Element segment;
  } vertex_properties_t;

  // Create a struct to hold properties for each edge
  typedef struct edge_properties {
    Element facet;
    Real weight;
  } edge_properties_t;

  // Define the type of the graph
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                                vertex_properties_t, edge_properties_t>
      graph_t;
  typedef typename graph_t::vertex_descriptor vertex_descriptor_t;
  typedef typename graph_t::edge_descriptor edge_descriptor_t;
  graph_t g;

  // prepare the list of facets sharing the central node
  CSR<Element> facets_to_nodes;
  MeshUtils::buildNode2Elements(mesh_facets, facets_to_nodes, dim - 1);
  auto && facets_to_node = facets_to_nodes.getRow(cent_node);
  std::set<Element> facets_added;
  std::map<Element, vertex_descriptor_t> segment_vertex_map;

  // add starting and ending segment to the graph
  vertex_descriptor_t v_start = boost::add_vertex(g);
  g[v_start].segment = limit_segments(0);
  vertex_descriptor_t v_end = boost::add_vertex(g);
  g[v_end].segment = limit_segments(1);
  segment_vertex_map[limit_segments(0)] = v_start;
  segment_vertex_map[limit_segments(1)] = v_end;

  for (auto && facet : facets_to_node) {
    // discard undesired facets first
    if (facet == limit_facets(0) or facet == limit_facets(1))
      continue;
    if (not isFacetAndNodesGood(facet))
      continue;
    if (belong2SameElement(limit_facets(0), facet) or
        belong2SameElement(limit_facets(1), facet))
      continue;

    // check the effective stress on this facet
    auto & facet_material_index =
        coh_model.getFacetMaterial(facet.type, facet.ghost_type)(facet.element);
    auto & mat = model.getMaterial(facet_material_index);
    auto * mat_coh = dynamic_cast<MaterialCohesive *>(&mat);
    if (mat_coh == nullptr)
      continue;
    auto && mat_facet_filter =
        mat_coh->getFacetFilter(facet.type, facet.ghost_type);
    UInt facet_local_num = mat_facet_filter.find(facet.element);
    AKANTU_DEBUG_ASSERT(facet_local_num != UInt(-1),
                        "Local facet number was not identified");
    if (not mat_coh->isInternal<Real>("effective_stresses", _ek_regular))
      continue;
    auto & eff_stress_array = mat.getInternal<Real>("effective_stresses")(
        facet.type, facet.ghost_type);
    // eff stress is always non-negative
    auto eff_stress = eff_stress_array(facet_local_num);
    if (eff_stress == 0)
      eff_stress = 1e-3;
    // if (eff_stress < 1)
    //   continue;

    // identify 2 sides of a facets
    Array<Element> side_segments;
    auto && facet_segments = mesh_facets.getSubelementToElement(facet);
    for (auto & segment : facet_segments) {
      auto ret = std::find(segments_to_point.begin(), segments_to_point.end(),
                           segment);
      if (ret != segments_to_point.end()) {
        side_segments.push_back(segment);
      }
    }
    AKANTU_DEBUG_ASSERT(side_segments.size() == 2,
                        "Number of side facets is wrong");

    // add corresponding vertex (if not there)and edge into graph
    Array<vertex_descriptor_t> v(2);
    for (UInt i = 0; i < 2; i++) {
      auto it = segment_vertex_map.find(side_segments(i));
      if (it != segment_vertex_map.end()) {
        v(i) = it->second;
      } else {
        v(i) = boost::add_vertex(g);
        g[v(i)].segment = side_segments(i);
        segment_vertex_map[side_segments(i)] = v(i);
      }
    }

    // add corresponding edge into the graph
    auto ret = facets_added.emplace(facet);
    if (ret.second) {
      // compute facet's area to use as a weight

      // add edge
      std::pair<edge_descriptor_t, bool> e = boost::add_edge(v(0), v(1), g);
      g[e.first].weight = 1 / eff_stress;
      g[e.first].facet = facet;
    }
  }

  // BGL Dijkstra's Shortest Paths here...
  std::vector<Real> distances(boost::num_vertices(g));
  std::vector<vertex_descriptor_t> predecessors(boost::num_vertices(g));

  boost::dijkstra_shortest_paths(
      g, v_start,
      boost::weight_map(boost::get(&edge_properties_t::weight, g))
          .distance_map(boost::make_iterator_property_map(
              distances.begin(), boost::get(boost::vertex_index, g)))
          .predecessor_map(boost::make_iterator_property_map(
              predecessors.begin(), boost::get(boost::vertex_index, g))));

  // Extract the shortest path from v_end to v_start (reversed order)
  typedef std::vector<edge_descriptor_t> path_t;
  path_t path;

  vertex_descriptor_t v = v_end;
  for (vertex_descriptor_t u = predecessors[v]; u != v;
       v = u, u = predecessors[v]) {
    std::pair<edge_descriptor_t, bool> edge_pair = boost::edge(u, v, g);
    path.push_back(edge_pair.first);
  }

  // compute average effective stress over the loop and discard if <1
  Array<Element> facet_loop_tmp;
  for (auto edge : path) {
    facet_loop_tmp.push_back(g[edge].facet);
  }
  auto av_eff_stress = averageEffectiveStressInMultipleFacets(facet_loop_tmp);
  if (av_eff_stress < 1)
    return facets_in_loop;

  // check path for smoothness
  bool recompute_path{false};
  do {
    recompute_path = false;
    if (path.size()) {
      for (UInt i = 0; i <= path.size(); i++) {
        // evaluate two neighboring facets
        Element facet1, facet2;
        if (i == 0) {
          facet1 = limit_facets(1);
          facet2 = g[path[i]].facet;
        } else if (i == path.size()) {
          facet1 = g[path[i - 1]].facet;
          facet2 = limit_facets(0);
        } else {
          facet1 = g[path[i - 1]].facet;
          facet2 = g[path[i]].facet;
        }
        // discard short distances
        auto facet1_indiam =
            MeshUtils::getInscribedCircleDiameter(model, facet1);
        auto facet2_indiam =
            MeshUtils::getInscribedCircleDiameter(model, facet2);
        auto hypotenuse = sqrt(facet1_indiam * facet1_indiam +
                               facet2_indiam * facet2_indiam) /
                          2;
        auto dist = MeshUtils::distanceBetweenTwoFacetsAlligned(mesh_facets,
                                                                facet1, facet2);
        // get abs of dot product between two normals
        Real dot =
            MeshUtils::cosSharpAngleBetween2Facets(model, facet1, facet2);

        if (dist < hypotenuse or dot < min_dot) {
          if (i < path.size()) {
            boost::remove_edge(path[i], g);
          } else {
            boost::remove_edge(path[i - 1], g);
          }
          recompute_path = true;
          break;
        }
      }
      // find new shortest path without problematic edge
      if (recompute_path) {
        boost::dijkstra_shortest_paths(
            g, v_start,
            boost::weight_map(boost::get(&edge_properties_t::weight, g))
                .distance_map(boost::make_iterator_property_map(
                    distances.begin(), boost::get(boost::vertex_index, g)))
                .predecessor_map(boost::make_iterator_property_map(
                    predecessors.begin(), boost::get(boost::vertex_index, g))));
        path.clear();
        v = v_end;
        for (vertex_descriptor_t u = predecessors[v]; u != v;
             v = u, u = predecessors[v]) {
          std::pair<edge_descriptor_t, bool> edge_pair = boost::edge(u, v, g);
          path.push_back(edge_pair.first);
        }
      }
    }
  } while (recompute_path);

  // add the shortest path facets
  for (typename path_t::reverse_iterator riter = path.rbegin();
       riter != path.rend(); ++riter) {
    vertex_descriptor_t u_tmp = boost::source(*riter, g);
    vertex_descriptor_t v_tmp = boost::target(*riter, g);
    edge_descriptor_t e_tmp = boost::edge(u_tmp, v_tmp, g).first;

    facets_in_loop.push_back(g[e_tmp].facet);
  }
  return facets_in_loop;
}
/* ------------------------------------------------------------------- */
Array<Element>
RVETools::findSingleFacetLoop(const Array<Element> & limit_facets,
                              const Array<Element> & limit_segments,
                              const UInt & cent_node) {

  auto & mesh = model.getMesh();
  auto & mesh_facets = mesh.getMeshFacets();
  auto dim = mesh.getSpatialDimension();

  AKANTU_DEBUG_ASSERT(limit_facets(0) != limit_facets(1),
                      "Limit facets are equal");
  AKANTU_DEBUG_ASSERT(limit_segments(0) != limit_segments(1),
                      "Limit segments are equal");
  AKANTU_DEBUG_ASSERT(
      dim == 3, "findFacetsLoopFromSegment2Segment currently supports only 3D");
  // get the point el corresponding to the central node
  CSR<Element> points_to_nodes;
  MeshUtils::buildNode2Elements(mesh_facets, points_to_nodes, 0);
  auto && points_to_node = points_to_nodes.getRow(cent_node);
  auto point_el = *points_to_node.begin();
  auto & segments_to_node = mesh_facets.getElementToSubelement(point_el);

  AKANTU_DEBUG_ASSERT(point_el.type == _point_1,
                      "Provided node element type is wrong");

  for (UInt i = 0; i < 2; i++) {
    // check the dimensions
    AKANTU_DEBUG_ASSERT(
        Mesh::getSpatialDimension(limit_facets(i).type) == dim - 1,
        "Facet" << limit_facets(i) << " has wrong spatial dimension");
    AKANTU_DEBUG_ASSERT(
        Mesh::getSpatialDimension(limit_segments(i).type) == dim - 2,
        "Segment" << limit_segments(i) << " has wrong spatial dimension");

    // check if the provided node is shared by two segments
    auto ret = std::find(segments_to_node.begin(), segments_to_node.end(),
                         limit_segments(i));
    AKANTU_DEBUG_ASSERT(ret != segments_to_node.end(),
                        "Segment " << limit_segments(i)
                                   << " doesn't contain the central node");

    // check if limit facets contain limit segments
    auto & facets_to_segment =
        mesh_facets.getElementToSubelement(limit_segments(i));
    auto ret1 = std::find(facets_to_segment.begin(), facets_to_segment.end(),
                          limit_facets(i));
    AKANTU_DEBUG_ASSERT(ret1 != facets_to_segment.end(),
                        "Facet " << limit_facets(i)
                                 << " doesn't contain the segment "
                                 << limit_segments(i));
  }

  Array<Element> facets_in_loop;

  // check if two segments are not connected by a single facet
  auto st_segment_neighbors =
      mesh_facets.getElementToSubelement(limit_segments(0));
  auto end_segment_neighbors =
      mesh_facets.getElementToSubelement(limit_segments(1));
  std::sort(st_segment_neighbors.begin(), st_segment_neighbors.end());
  std::sort(end_segment_neighbors.begin(), end_segment_neighbors.end());

  Array<Element> shared_facets(st_segment_neighbors.size());
  // An iterator to the end of the constructed range.
  auto it = std::set_intersection(
      st_segment_neighbors.begin(), st_segment_neighbors.end(),
      end_segment_neighbors.begin(), end_segment_neighbors.end(),
      shared_facets.begin());
  shared_facets.resize(it - shared_facets.begin());

  // if there is a single facet connecting two segments test it
  if (shared_facets.size() == 1) {
    bool bad_facet{false};
    if (not isFacetAndNodesGood(shared_facets(0)))
      bad_facet = true;
    // check angle with starting and ending -> not to be sharp
    bool sharp_angle{false};
    auto facet_indiam =
        MeshUtils::getInscribedCircleDiameter(model, shared_facets(0));
    for (auto && limit_facet : limit_facets) {
      auto limit_facet_indiam =
          MeshUtils::getInscribedCircleDiameter(model, limit_facet);
      auto hypotenuse = sqrt(limit_facet_indiam * limit_facet_indiam +
                             facet_indiam * facet_indiam) /
                        2;
      auto dist = MeshUtils::distanceBetweenTwoFacetsAlligned(
          mesh_facets, limit_facet, shared_facets(0));

      if (dist <= hypotenuse) {
        sharp_angle = true;
        break;
      }
    }
    if (not sharp_angle and not bad_facet) {
      facets_in_loop.push_back(shared_facets(0));
    }
  }
  return facets_in_loop;
}
/* ------------------------------------------------------------------- */
bool RVETools::isFacetAndNodesGood(const Element & facet, UInt material_id,
                                   bool check_crack_nodes) {

  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  auto & inserter = coh_model.getElementInserter();
  auto facet_type = facet.type;
  auto gt = facet.ghost_type;
  auto & mesh_facets = inserter.getMeshFacets();
  Vector<UInt> facet_nodes = mesh_facets.getConnectivity(facet);
  auto & check_facets = inserter.getCheckFacets(facet_type, gt);

  // check if the facet is not ghost
  if (gt == _ghost)
    return false;
  // check if the facet is allowed to be inserted
  if (check_facets(facet.element) == false)
    return false;
  for (auto & facet_node : facet_nodes) {
    // check if the facet's nodes are not pure ghost nodes
    if (model.getMesh().isPureGhostNode(facet_node))
      return false;

    // check if the facet's nodes are the crack nodes
    if (check_crack_nodes and this->crack_nodes(facet_node))
      return false;
  }

  // check if the facet material is good
  if (material_id != UInt(-1)) {
    const auto & facet_material =
        coh_model.getFacetMaterial(facet_type, gt)(facet.element);
    if (facet_material != material_id)
      return false;
  }

  return true;
}
/* ------------------------------------------------------------------- */
bool RVETools::belong2SameElement(const Element & facet1,
                                  const Element & facet2) {

  auto & mesh = model.getMesh();
  auto & mesh_facets = mesh.getMeshFacets();

  Array<Element> f1_neighbors, f2_neighbors;
  auto facet1_neighbors = mesh_facets.getElementToSubelement(facet1);
  for (auto && neighbor : facet1_neighbors) {
    if (mesh.getKind(neighbor.type) == _ek_regular)
      f1_neighbors.push_back(neighbor);
    if (mesh.getKind(neighbor.type) == _ek_cohesive) {
      auto && coh_neighbors = mesh_facets.getSubelementToElement(neighbor);
      for (auto && coh_neighbor : coh_neighbors) {
        if (coh_neighbor == facet1)
          continue;
        auto coh_facet_neighbors =
            mesh_facets.getElementToSubelement(coh_neighbor);
        for (auto && coh_facet_neighbor : coh_facet_neighbors) {
          if (mesh.getKind(coh_facet_neighbor.type) == _ek_regular)
            f1_neighbors.push_back(coh_facet_neighbor);
        }
      }
    }
  }

  auto facet2_neighbors = mesh_facets.getElementToSubelement(facet2);
  for (auto && neighbor : facet2_neighbors) {
    if (mesh.getKind(neighbor.type) == _ek_regular)
      f2_neighbors.push_back(neighbor);
    if (mesh.getKind(neighbor.type) == _ek_cohesive) {
      auto && coh_neighbors = mesh_facets.getSubelementToElement(neighbor);
      for (auto && coh_neighbor : coh_neighbors) {
        if (coh_neighbor == facet2)
          continue;
        auto coh_facet_neighbors =
            mesh_facets.getElementToSubelement(coh_neighbor);
        for (auto && coh_facet_neighbor : coh_facet_neighbors) {
          if (mesh.getKind(coh_facet_neighbor.type) == _ek_regular)
            f2_neighbors.push_back(coh_facet_neighbor);
        }
      }
    }
  }
  std::sort(f1_neighbors.begin(), f1_neighbors.end());
  std::sort(f2_neighbors.begin(), f2_neighbors.end());

  Array<Element> shared_neighbors(f1_neighbors.size());
  // An iterator to the end of the constructed range.
  auto it = std::set_intersection(f1_neighbors.begin(), f1_neighbors.end(),
                                  f2_neighbors.begin(), f2_neighbors.end(),
                                  shared_neighbors.begin());
  shared_neighbors.resize(it - shared_neighbors.begin());

  if (shared_neighbors.size() == 0)
    return false;

  // check kind of shared neighbor
  auto shared_neighbor = shared_neighbors(0);
  if (mesh.getKind(shared_neighbor.type) != _ek_regular)
    return false;

  return true;
}
/* ------------------------------------------------------------------- */
bool RVETools::isNodeWithinMaterial(Element & node, UInt material_id) {

  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  auto & inserter = coh_model.getElementInserter();
  auto & mesh_facets = inserter.getMeshFacets();
  auto dim = coh_model.getSpatialDimension();

  CSR<Element> nodes_to_facets;
  MeshUtils::buildNode2Elements(mesh_facets, nodes_to_facets, dim - 1);

  for (auto & facet : nodes_to_facets.getRow(node.element)) {
    auto & facet_material =
        coh_model.getFacetMaterial(facet.type, facet.ghost_type)(facet.element);
    if (facet_material != material_id)
      return false;
  }

  return true;
}
/* ----------------------------------------------------------------------- */
bool RVETools::pickFacetNeighbors(Element & cent_facet) {
  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  auto & inserter = coh_model.getElementInserter();
  auto & mesh = model.getMesh();
  auto & mesh_facets = inserter.getMeshFacets();
  auto dim = mesh.getSpatialDimension();
  auto facet_type = cent_facet.type;
  auto facet_gt = cent_facet.ghost_type;
  auto facet_nb = cent_facet.element;
  auto & doubled_facets = mesh_facets.getElementGroup("doubled_facets");
  auto & doubled_nodes = mesh.getNodeGroup("doubled_nodes");
  auto & facet_conn = mesh_facets.getConnectivity(facet_type, facet_gt);
  auto & facet_material_index =
      coh_model.getFacetMaterial(facet_type, facet_gt)(facet_nb);

  // get list of the subel connected to the facet (nodes in 2D, segments
  // in 3D)
  const Vector<Element> & subelements_to_element =
      mesh_facets.getSubelementToElement(cent_facet);

  Array<Element> neighbors(subelements_to_element.size());
  neighbors.set(ElementNull);
  // max weight parameter computed for each subelement iteration
  // here this parameter = dot product between two normals
  // minimum 60 deg between normals -> 120 between dips
  Real weight_parameter_threshold(0.5);
  std::vector<Real> max_weight_parameters(subelements_to_element.size(),
                                          weight_parameter_threshold);

  // identify a facet's neighbor through each subelement
  for (UInt i : arange(subelements_to_element.size())) {
    auto & subel = subelements_to_element(i);
    auto & connected_elements = mesh_facets.getElementToSubelement(
        subel.type, subel.ghost_type)(subel.element);

    for (auto & connected_element : connected_elements) {
      // check all possible unsatisfactory conditions
      if (connected_element.type != facet_type)
        continue;
      if (connected_element.element == facet_nb)
        continue;
      auto & candidate_facet_material_index = coh_model.getFacetMaterial(
          connected_element.type,
          connected_element.ghost_type)(connected_element.element);
      if (candidate_facet_material_index != facet_material_index)
        continue;
      if (not inserter.getCheckFacets(connected_element.type,
                                      connected_element.ghost_type)(
              connected_element.element))
        continue;
      // discard facets intersecting other crack elements
      bool crack_node{false};
      for (UInt j : arange(facet_conn.getNbComponent())) {
        auto node = facet_conn(connected_element.element, j);
        if (this->crack_nodes(node))
          crack_node = true;
      }
      if (crack_node)
        continue;

      // get inscribed diameter
      Real facet_indiam =
          MeshUtils::getInscribedCircleDiameter(model, cent_facet);

      // get distance between two barycenters
      auto dist = MeshUtils::distanceBetweenTwoFacetsAlligned(
          mesh_facets, cent_facet, connected_element);

      // ad-hoc rule on barycenters spacing
      // it should discard all elements under sharp angle
      if (dist < facet_indiam)
        continue;

      // get abs of dot product between two normals
      Real dot = MeshUtils::cosSharpAngleBetween2Facets(model, cent_facet,
                                                        connected_element);

      auto weight_parameter = dot;
      if (weight_parameter > max_weight_parameters[i]) {
        max_weight_parameters[i] = weight_parameter;
        neighbors(i) = connected_element;
      }
    }
  }

  // different insertion procedures for 2D and 3D cases
  switch (dim) {
  case 2: {
    if (neighbors.find(ElementNull) == UInt(-1)) {
      for (auto & neighbor : neighbors) {
        doubled_facets.add(neighbor);
        for (UInt node : arange(2)) {
          this->crack_nodes(facet_conn(neighbor.element, node)) = true;
        }
      }
      return true;
    } else
      return false;
  }
  case 3: {
    auto max_el_pos = std::max_element(max_weight_parameters.begin(),
                                       max_weight_parameters.end());
    if (*max_el_pos > weight_parameter_threshold) {
      auto max_param_index =
          std::distance(max_weight_parameters.begin(), max_el_pos);
      auto & neighbor = neighbors(max_param_index);
      doubled_facets.add(neighbor);
      for (UInt node : arange(facet_conn.getNbComponent())) {
        doubled_nodes.add(facet_conn(neighbor.element, node));
        this->crack_nodes(facet_conn(neighbor.element, node)) = true;
      }
      return true;
    } else
      return false;
  }
  default:
    AKANTU_EXCEPTION("Provided dimension is not supported");
  }
}
/* ------------------------------------------------------------------ */
void RVETools::preventCohesiveInsertionInNeighbors() {
  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  auto & inserter = coh_model.getElementInserter();
  auto & mesh = model.getMesh();
  auto & mesh_facets = inserter.getMeshFacets();
  auto dim = mesh.getSpatialDimension();
  const GhostType gt = _not_ghost;
  auto el_type = *mesh.elementTypes(dim, gt, _ek_regular).begin();
  auto facet_type = Mesh::getFacetType(el_type);

  auto & subfacets_to_facets =
      mesh_facets.getSubelementToElement(facet_type, gt);
  auto nb_subfacet = subfacets_to_facets.getNbComponent();
  auto subf_to_fac_it = subfacets_to_facets.begin(nb_subfacet);

  auto & el_group = mesh_facets.getElementGroup("doubled_facets");
  Array<UInt> element_ids = el_group.getElements(facet_type);
  for (UInt i : arange(element_ids.size() / 2)) {
    auto facet1 = element_ids(i);
    auto facet2 = element_ids(i + element_ids.size() / 2);
    std::vector<Element> subfacets1(nb_subfacet);
    std::vector<Element> subfacets2(nb_subfacet);
    for (UInt i : arange(nb_subfacet)) {
      subfacets1[i] = subf_to_fac_it[facet1](i);
      subfacets2[i] = subf_to_fac_it[facet2](i);
    }
    std::sort(subfacets1.begin(), subfacets1.end());
    std::sort(subfacets2.begin(), subfacets2.end());
    std::vector<Element> dif_subfacets(2 * nb_subfacet);
    std::vector<Element>::iterator it;
    // identify subfacets not shared by two facets
    it = std::set_difference(subfacets1.begin(), subfacets1.end(),
                             subfacets2.begin(), subfacets2.end(),
                             dif_subfacets.begin());
    dif_subfacets.resize(it - dif_subfacets.begin());

    // prevent insertion of cohesives in the facets including these subf
    for (auto & subfacet : dif_subfacets) {
      auto neighbor_facets = mesh_facets.getElementToSubelement(subfacet);
      for (auto & neighbor_facet : neighbor_facets) {
        inserter.getCheckFacets(facet_type, gt)(neighbor_facet.element) = false;
      }
    }
  }
}

/* ------------------------------------------------------------------- */
void RVETools::insertOppositeFacetsAndCohesives() {
  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  auto & inserter = coh_model.getElementInserter();
  auto & insertion = inserter.getInsertionFacetsByElement();
  auto & mesh = model.getMesh();
  auto & mesh_facets = inserter.getMeshFacets();
  auto dim = mesh.getSpatialDimension();
  const GhostType gt = _not_ghost;
  const auto & pos = mesh.getNodes();
  const auto pos_it = make_view(pos, dim).begin();
  auto & doubled_facets = mesh_facets.getElementGroup("doubled_facets");

  // sort facet numbers in doubled_facets_group to comply with new_elements
  // event passed by the inserter
  doubled_facets.optimize();

  // instruct the inserter which facets to duplicate
  for (auto & type : doubled_facets.elementTypes(dim - 1)) {
    const auto element_ids = doubled_facets.getElements(type, gt);
    // iterate over facets to duplicate
    for (auto && el : element_ids) {
      inserter.getCheckFacets(type, gt)(el) = false;
      Element new_facet{type, el, gt};
      insertion(new_facet) = true;
    }
  }

  // duplicate facets and insert coh els
  inserter.insertElements(false);

  // add facets connectivity to the mesh and the element group
  NewElementsEvent new_facets_event;
  for (auto & type : doubled_facets.elementTypes(dim - 1)) {
    const auto element_ids = doubled_facets.getElements(type, gt);

    if (not mesh.getConnectivities().exists(type, gt))
      mesh.addConnectivityType(type, gt);
    auto & facet_conn = mesh.getConnectivity(type, gt);
    auto & mesh_facet_conn = mesh_facets.getConnectivity(type, gt);
    const UInt nb_nodes_facet = facet_conn.getNbComponent();
    const auto mesh_facet_conn_it =
        make_view(mesh_facet_conn, nb_nodes_facet).begin();

    // iterate over duplicated facets
    for (auto && el : element_ids) {
      Vector<UInt> facet_nodes = mesh_facet_conn_it[el];
      facet_conn.push_back(facet_nodes);
      Element new_facet{type, facet_conn.size() - 1, gt};
      new_facets_event.getList().push_back(new_facet);
    }
  }

  MeshUtils::fillElementToSubElementsData(mesh);
  mesh.sendEvent(new_facets_event);

  // fill up crack facets vector
  this->crack_facets_from_mesh.resize(this->crack_facets.size());
  for (UInt i = 0; i != this->crack_facets.size(); i++) {
    auto single_site_facets = this->crack_facets(i);
    for (auto & facet : single_site_facets) {
      Element cohesive{ElementNull};
      auto & connected_els = mesh_facets.getElementToSubelement(facet);
      for (auto & connected_el : connected_els) {
        if (mesh.getKind(connected_el.type) == _ek_cohesive) {
          cohesive = connected_el;
          break;
        }
      }
      AKANTU_DEBUG_ASSERT(cohesive != ElementNull,
                          "Connected cohesive element not identified");
      auto && connected_subels = mesh.getSubelementToElement(cohesive);
      for (auto & connected_subel : connected_subels) {
        if (connected_subel.type != facet.type)
          continue;
        this->crack_facets_from_mesh(i).push_back(connected_subel);
      }
    }
  }

  // create an element group with nodes to apply Dirichlet
  model.getMesh().createElementGroupFromNodeGroup("doubled_nodes",
                                                  "doubled_nodes", dim - 1);

  // update FEEngineBoundary with new elements
  model.getFEEngineBoundary().initShapeFunctions(_not_ghost);
  model.getFEEngineBoundary().initShapeFunctions(_ghost);
  model.getFEEngineBoundary().computeNormalsOnIntegrationPoints(_not_ghost);
  model.getFEEngineBoundary().computeNormalsOnIntegrationPoints(_ghost);
}
/* ----------------------------------------------------------------------- */
void RVETools::assignCrackNumbers() {
  auto & mesh = model.getMesh();
  auto & mesh_facets = mesh.getMeshFacets();
  auto dim = mesh.getSpatialDimension();
  const GhostType gt = _not_ghost;

  for (auto && type_facet : mesh_facets.elementTypes(dim - 1)) {
    ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);
    auto & crack_numbers =
        mesh.getDataPointer<UInt>("crack_numbers", type_cohesive);
    auto & coh_conn = mesh.getConnectivity(type_cohesive, gt);
    auto nb_coh_nodes = coh_conn.getNbComponent();
    auto coh_conn_it = coh_conn.begin(coh_conn.getNbComponent());

    // initialize a graph
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>
        Graph;
    Graph graph;

    // build the graph based on the connectivity of cohesive elements
    for (UInt i = 0; i < coh_conn.size(); i++) {
      for (UInt j = i + 1; j < coh_conn.size(); j++) {
        auto nodes_1 = coh_conn_it[i];
        auto nodes_2 = coh_conn_it[j];
        std::vector<UInt> nod_1(nb_coh_nodes);
        std::vector<UInt> nod_2(nb_coh_nodes);
        for (UInt k : arange(nb_coh_nodes)) {
          nod_1[k] = nodes_1(k);
          nod_2[k] = nodes_2(k);
        }

        std::sort(nod_1.begin(), nod_1.end());
        std::sort(nod_2.begin(), nod_2.end());
        std::vector<UInt> common_nodes(nb_coh_nodes);
        std::vector<UInt>::iterator it;

        // An iterator to the end of the constructed range.
        it = std::set_intersection(nod_1.begin(), nod_1.end(), nod_2.begin(),
                                   nod_2.end(), common_nodes.begin());
        common_nodes.resize(it - common_nodes.begin());

        if (common_nodes.size() > 1) {
          boost::add_edge(i, j, graph);
        }
      }
    }

    // connectivity and number of components of a graph
    std::vector<int> component(boost::num_vertices(graph));
    int num = boost::connected_components(graph, &component[0]);

    // shift the crack numbers by the number of cracks on processors with the
    // lower rank
    auto && comm = akantu::Communicator::getWorldCommunicator();
    comm.exclusiveScan(num);
    for (auto & component_nb : component) {
      component_nb += num;
    }

    // assign a corresponding crack flag
    for (UInt coh_el : arange(component.size())) {
      crack_numbers(coh_el) = component[coh_el];
    }
  }
}
/* ------------------------------------------------------------------- */
void RVETools::insertGap(const Real gap_ratio) {
  AKANTU_DEBUG_IN();
  if (gap_ratio == 0)
    return;

  auto & mesh = model.getMesh();
  MeshAccessor mesh_accessor(mesh);
  auto & pos2modify = mesh_accessor.getNodes();

  auto dim = mesh.getSpatialDimension();
  const GhostType gt = _not_ghost;
  auto & el_group = mesh.getElementGroup("doubled_nodes");
  auto & pos = mesh.getNodes();
  auto pos_it = make_view(pos, dim).begin();
  auto pos2modify_it = make_view(pos2modify, dim).begin();

  for (auto & type : el_group.elementTypes(dim - 1)) {
    auto & facet_conn = mesh.getConnectivity(type, gt);
    const UInt nb_nodes_facet = facet_conn.getNbComponent();
    auto facet_nodes_it = make_view(facet_conn, nb_nodes_facet).begin();
    const auto element_ids = el_group.getElements(type, gt);
    auto && fe_engine = model.getFEEngineBoundary();
    auto nb_qpoints_per_facet = fe_engine.getNbIntegrationPoints(type, gt);
    const auto & normals_on_quad =
        fe_engine.getNormalsOnIntegrationPoints(type, gt);
    auto normals_it = make_view(normals_on_quad, dim).begin();
    for (UInt i : arange(element_ids.size())) {
      auto el_id = element_ids(i);
      Vector<Real> normal = normals_it[el_id * nb_qpoints_per_facet];
      if (i < element_ids.size() / 2)
        normal *= -1;
      // compute segment length
      auto facet_nodes = facet_nodes_it[el_id];
      UInt A, B;
      A = facet_nodes(0);
      B = facet_nodes(1);
      Vector<Real> AB = Vector<Real>(pos_it[B]) - Vector<Real>(pos_it[A]);
      Real correction = AB.norm() * gap_ratio / 2;
      Vector<Real> half_opening_vector = normal * correction;
      for (auto j : arange(nb_nodes_facet)) {
        Vector<Real> node_pos(pos2modify_it[facet_nodes(j)]);
        node_pos += half_opening_vector;
        this->modified_pos(facet_nodes(j)) = true;
      }
    }
  }

  // update FEEngine & FEEngineBoundary with new elements
  model.getFEEngine().initShapeFunctions(_not_ghost);
  model.getFEEngine().initShapeFunctions(_ghost);
  model.getFEEngineBoundary().initShapeFunctions(_not_ghost);
  model.getFEEngineBoundary().initShapeFunctions(_ghost);
  model.getFEEngineBoundary().computeNormalsOnIntegrationPoints(_not_ghost);
  model.getFEEngineBoundary().computeNormalsOnIntegrationPoints(_ghost);

  AKANTU_DEBUG_OUT();
}
/* ------------------------------------------------------------------- */
void RVETools::insertGap3D(const Real gap_ratio) {
  AKANTU_DEBUG_IN();
  if (gap_ratio == 0)
    return;

  auto & mesh = model.getMesh();
  auto & mesh_facets = mesh.getMeshFacets();
  MeshAccessor mesh_accessor(mesh);
  auto & pos2modify = mesh_accessor.getNodes();

  auto dim = mesh.getSpatialDimension();
  const GhostType gt = _not_ghost;
  auto & el_group = mesh_facets.getElementGroup("doubled_facets");
  auto & pos = mesh.getNodes();
  auto pos_it = make_view(pos, dim).begin();
  auto pos2modify_it = make_view(pos2modify, dim).begin();

  for (auto & type : el_group.elementTypes(dim - 1)) {
    if ((dim == 3) and (type != _triangle_6))
      AKANTU_EXCEPTION("The only facet type supported in 3D is _triangle_6");
    const Array<UInt> element_ids = el_group.getElements(type, gt);
    auto && fe_engine = model.getFEEngineBoundary();
    auto nb_qpoints_per_facet = fe_engine.getNbIntegrationPoints(type, gt);

    // identify pairs of facets
    Array<UInt> facets_considered;
    for (UInt i : arange(element_ids.size())) {
      auto el_id = element_ids(i);
      Element element{type, el_id, gt};

      // if element was already considered -> skip
      if (facets_considered.find(el_id) != UInt(-1))
        continue;

      Array<UInt> element_ids_half(element_ids.size() / 2);
      UInt starting_pos;
      if (i < element_ids.size() / 2) {
        starting_pos = 0;
      } else {
        starting_pos = element_ids.size() / 2;
      }
      for (auto j : arange(element_ids_half.size())) {
        element_ids_half(j) = element_ids(j + starting_pos);
      }

      // find a neighbor
      const Vector<Element> & subelements_to_element =
          mesh_facets.getSubelementToElement(element);
      // identify a facet's neighbor through each subelement
      Element neighbor(ElementNull);
      Element border(ElementNull);
      for (UInt i : arange(subelements_to_element.size())) {
        auto & subel = subelements_to_element(i);
        auto & connected_elements = mesh_facets.getElementToSubelement(subel);
        for (auto & connected_element : connected_elements) {
          // check all possible unsatisfactory conditions
          if (connected_element.type != type)
            continue;
          if (connected_element.element == el_id)
            continue;
          // search for this neighbor in the doubled_nodes el group
          UInt neighbor_pos_in_el_ids =
              element_ids_half.find(connected_element.element);
          if (neighbor_pos_in_el_ids == UInt(-1))
            continue;

          // if found -> subelement is the border
          neighbor = connected_element;
          border = subel;
          facets_considered.push_back(el_id);
          facets_considered.push_back(connected_element.element);
          break;
        }
        if (neighbor != ElementNull)
          break;
      }
      AKANTU_DEBUG_ASSERT(neighbor != ElementNull,
                          "Neighbor for the facet was not identified for the "
                          "purpose of the gap insertion");

      // average the normal in between two neighbors
      const auto & normals_on_quad =
          fe_engine.getNormalsOnIntegrationPoints(type, gt);
      auto normals_it = make_view(normals_on_quad, dim).begin();
      Vector<Real> normal = normals_it[el_id * nb_qpoints_per_facet];
      normal +=
          Vector<Real>(normals_it[neighbor.element * nb_qpoints_per_facet]);
      normal /= normal.norm();
      // and correct it for the direction
      if (i < element_ids.size() / 2) {
        normal *= -1;
      }

      // compute segment length
      AKANTU_DEBUG_ASSERT(
          border.type == _segment_3,
          "The only supported segment type in 3D is _segment_3");
      auto & segment_conn =
          mesh_facets.getConnectivity(border.type, border.ghost_type);
      const UInt nb_nodes_segment = segment_conn.getNbComponent();
      auto segment_nodes_it = make_view(segment_conn, nb_nodes_segment).begin();
      auto segment_nodes = segment_nodes_it[border.element];
      // 2 end nodes and the middle node
      UInt A, B, middle_node;
      A = segment_nodes(0);
      B = segment_nodes(1);
      middle_node = segment_nodes(2);
      Vector<Real> AB = Vector<Real>(pos_it[B]) - Vector<Real>(pos_it[A]);
      Real correction = AB.norm() * gap_ratio / 2;
      Vector<Real> half_opening_vector = normal * correction;
      Vector<Real> node_pos(pos2modify_it[middle_node]);
      node_pos += half_opening_vector;
      this->modified_pos(middle_node) = true;
    }
  }

  // update FEEngine & FEEngineBoundary with new elements
  model.getFEEngine().initShapeFunctions(_not_ghost);
  model.getFEEngine().initShapeFunctions(_ghost);
  model.getFEEngineBoundary().initShapeFunctions(_not_ghost);
  model.getFEEngineBoundary().initShapeFunctions(_ghost);
  model.getFEEngineBoundary().computeNormalsOnIntegrationPoints(_not_ghost);
  model.getFEEngineBoundary().computeNormalsOnIntegrationPoints(_ghost);

  AKANTU_DEBUG_OUT();
}

/* ------------------------------------------------------------------- */
template <UInt dim>
UInt RVETools::insertLoopsOfStressedCohesives(Real min_dot) {
  AKANTU_DEBUG_IN();

  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  CohesiveElementInserter & inserter = coh_model.getElementInserter();

  if (not coh_model.getIsExtrinsic()) {
    AKANTU_EXCEPTION(
        "This function can only be used for extrinsic cohesive elements");
  }

  coh_model.interpolateStress();

  // compute effective stress in all cohesive material linear sequent.
  for (auto && mat : model.getMaterials()) {
    auto * mat_coh =
        dynamic_cast<MaterialCohesiveLinearSequential<dim> *>(&mat);

    if (mat_coh == nullptr)
      continue;

    mat_coh->computeEffectiveStresses();
  }

  auto && mesh = coh_model.getMesh();
  const Mesh & mesh_facets = coh_model.getMeshFacets();
  UInt nb_new_elements(0);
  auto type_facet = *mesh_facets.elementTypes(dim - 1).begin();
  ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);
  auto & crack_numbers = mesh.getData<UInt>("crack_numbers", type_cohesive);

  // find global crack topology
  auto data = determineCrackSurface();
  auto & contour_subfacets_coh_el = std::get<0>(data);
  auto & surface_subfacets_crack_nb = std::get<1>(data);
  auto & surface_nodes = std::get<2>(data);
  auto & contour_nodes_subfacets = std::get<3>(data);

  // identify facets on contour where cohesives should be inserted
  std::map<UInt, std::map<UInt, UInt>> mat_index_facet_nbs_crack_nbs =
      findStressedFacetsLoops(contour_subfacets_coh_el,
                              surface_subfacets_crack_nb,
                              contour_nodes_subfacets, surface_nodes, min_dot);

  // insert cohesives depending on a material
  std::map<UInt, UInt> facet_nbs_crack_nbs_global;
  for (auto & pair : mat_index_facet_nbs_crack_nbs) {
    auto mat_index = pair.first;
    auto & facet_nbs_crack_nbs = pair.second;
    auto & mat = model.getMaterial(mat_index);
    auto * mat_coh =
        dynamic_cast<MaterialCohesiveLinearSequential<dim> *>(&mat);

    if (mat_coh == nullptr)
      continue;

    mat_coh->insertCohesiveElements(facet_nbs_crack_nbs, type_facet, false);

    facet_nbs_crack_nbs_global.insert(facet_nbs_crack_nbs.begin(),
                                      facet_nbs_crack_nbs.end());
  }

  // extend crack numbers array
  for (auto && pair : facet_nbs_crack_nbs_global) {
    auto && crack_nb = pair.second;
    crack_numbers.push_back(crack_nb);
  }

  // insert the most stressed cohesive element
  nb_new_elements = inserter.insertElements();

  // communicate crack numbers
  communicateCrackNumbers();

  AKANTU_DEBUG_OUT();
  return nb_new_elements;
}

template UInt RVETools::insertLoopsOfStressedCohesives<2>(Real min_dot);
template UInt RVETools::insertLoopsOfStressedCohesives<3>(Real min_dot);

/* ------------------------------------------------------------------- */
template <UInt dim> UInt RVETools::insertStressedCohesives() {
  AKANTU_DEBUG_IN();

  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  CohesiveElementInserter & inserter = coh_model.getElementInserter();

  if (not coh_model.getIsExtrinsic()) {
    AKANTU_EXCEPTION(
        "This function can only be used for extrinsic cohesive elements");
  }

  coh_model.interpolateStress();

  auto && mesh = coh_model.getMesh();
  const Mesh & mesh_facets = coh_model.getMeshFacets();
  UInt nb_new_elements(0);
  auto type_facet = *mesh_facets.elementTypes(dim - 1).begin();
  ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);
  auto & check_facets = inserter.getCheckFacets(type_facet);
  auto & crack_numbers = mesh.getData<UInt>("crack_numbers", type_cohesive);
  std::map<UInt, UInt> facet_nbs_crack_nbs;
  std::set<UInt> facet_nbs;

  // compute effective stress in all cohesive material linear sequent.
  for (auto && mat : model.getMaterials()) {
    auto * mat_coh =
        dynamic_cast<MaterialCohesiveLinearSequential<dim> *>(&mat);

    if (mat_coh == nullptr)
      continue;

    mat_coh->computeEffectiveStresses();

    // loop on each facet of the material
    auto && mat_facet_filter = mat_coh->getFacetFilter(type_facet);
    for (auto && data : enumerate(mat_facet_filter)) {
      auto facet_local_num = std::get<0>(data);
      auto facet_global_num = std::get<1>(data);
      if (check_facets(facet_global_num) == false)
        continue;

      auto & eff_stress_array =
          mat.getInternal<Real>("effective_stresses")(type_facet);

      // eff stress is always non-negative
      auto eff_stress = eff_stress_array(facet_local_num);
      if (eff_stress >= 1) {
        auto ret = facet_nbs.emplace(facet_global_num);
        if (ret.second) {
          facet_nbs_crack_nbs[facet_global_num] = 0;
        }
      }
    }
    mat_coh->insertCohesiveElements(facet_nbs_crack_nbs, type_facet, false);
  }

  // extend crack numbers array
  for (auto && pair : facet_nbs_crack_nbs) {
    auto && crack_nb = pair.second;
    crack_numbers.push_back(crack_nb);
  }

  // insert the most stressed cohesive element
  nb_new_elements = inserter.insertElements();

  // communicate crack numbers
  communicateCrackNumbers();

  AKANTU_DEBUG_OUT();
  return nb_new_elements;
}

template UInt RVETools::insertStressedCohesives<2>();
template UInt RVETools::insertStressedCohesives<3>();

/* --------------------------------------------------------------- */
std::map<UInt, std::map<UInt, UInt>> RVETools::findStressedFacetsLoops(
    const std::map<Element, Element> & contour_subfacets_coh_el,
    const std::map<Element, UInt> & /*surface_subfacets_crack_nb*/,
    const std::map<UInt, std::set<Element>> & contour_nodes_subfacets,
    const std::set<UInt> & /*surface_nodes*/, Real min_dot) {
  AKANTU_DEBUG_IN();

  const Mesh & mesh = model.getMesh();
  const Mesh & mesh_facets = mesh.getMeshFacets();
  auto dim = mesh.getSpatialDimension();
  std::map<UInt, std::map<UInt, UInt>> mat_index_facet_nbs_crack_nbs;
  std::map<UInt, std::pair<UInt, Array<Element>>> node_to_crack_nb_facet_loop;

  // clear the nodes_eff_stress array
  this->nodes_eff_stress.zero();

  // TODO : make iteration on multiple facet types
  auto type_facet = *mesh_facets.elementTypes(dim - 1).begin();
  auto nb_facets = mesh.getNbElement(type_facet);
  // skip if no facets of this type are present
  if (nb_facets) {
    for (auto && pair : contour_nodes_subfacets) {
      auto cent_node = pair.first;
      auto subfacets = pair.second;
      if (subfacets.size() != 2) {
        std::cout << "Number of contour segments connected to the node "
                  << cent_node << " is not equal to 2. Skipping" << std::endl;
        continue;
      }
      auto it = subfacets.begin();
      auto starting_segment(*it);
      std::advance(it, 1);
      auto ending_segment(*it);
      auto starting_coh_el = contour_subfacets_coh_el.at(starting_segment);
      auto ending_coh_el = contour_subfacets_coh_el.at(ending_segment);
      auto & crack_numbers = mesh.getData<UInt>(
          "crack_numbers", starting_coh_el.type, starting_coh_el.ghost_type);
      auto crack_nb = crack_numbers(starting_coh_el.element);

      auto starting_facet =
          mesh_facets.getSubelementToElement(starting_coh_el)(0);
      auto ending_facet = mesh_facets.getSubelementToElement(ending_coh_el)(0);

      Array<Element> limit_facets, limit_segments;
      limit_facets.push_back(starting_facet);
      limit_facets.push_back(ending_facet);
      limit_segments.push_back(starting_segment);
      limit_segments.push_back(ending_segment);

      auto facet_loop = findStressedFacetLoopAroundNode(
          limit_facets, limit_segments, cent_node, min_dot);
      // store the facet loop in the map
      node_to_crack_nb_facet_loop[cent_node] =
          std::make_pair(crack_nb, facet_loop);

      // evaluate average effective stress over facets of the loop
      auto av_eff_stress = averageEffectiveStressInMultipleFacets(facet_loop);

      // assign stress value to the node
      this->nodes_eff_stress(cent_node) = av_eff_stress;
    }
  }

  // communicate effective stresses by erasing the lowest values
  communicateEffStressesOnNodes();

  for (auto && pair : node_to_crack_nb_facet_loop) {
    auto && cent_node = pair.first;
    auto && crack_nb = pair.second.first;
    auto && facet_loop = pair.second.second;

    if (this->nodes_eff_stress(cent_node) > 0) {
      decomposeFacetLoopPerMaterial(facet_loop, crack_nb,
                                    mat_index_facet_nbs_crack_nbs);
    }
  }
  return mat_index_facet_nbs_crack_nbs;
}
/* --------------------------------------------------------------- */
Real RVETools::averageEffectiveStressInMultipleFacets(
    Array<Element> & facet_loop) {

  if (not facet_loop.size())
    return 0;

  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  const FEEngine & fe_engine = model.getFEEngine("FacetsFEEngine");

  Array<Real> eff_stress_array;
  ElementType type_facet;
  GhostType gt;
  Array<UInt> facet_nbs;

  for (auto & facet : facet_loop) {
    type_facet = facet.type;
    gt = facet.ghost_type;
    facet_nbs.push_back(facet.element);
    auto & facet_material_index =
        coh_model.getFacetMaterial(type_facet, gt)(facet.element);
    auto & mat = model.getMaterial(facet_material_index);
    auto * mat_coh = dynamic_cast<MaterialCohesive *>(&mat);

    if (mat_coh == nullptr)
      return 0;

    auto && mat_facet_filter = mat_coh->getFacetFilter(type_facet, gt);
    UInt facet_local_num = mat_facet_filter.find(facet.element);
    AKANTU_DEBUG_ASSERT(facet_local_num != UInt(-1),
                        "Local facet number was not identified");

    auto & eff_stress = mat.getInternal<Real>("effective_stresses")(
        type_facet, gt)(facet_local_num);
    for (__attribute__((unused)) UInt quad :
         arange(fe_engine.getNbIntegrationPoints(type_facet)))
      eff_stress_array.push_back(eff_stress);
  }

  Real int_eff_stress =
      fe_engine.integrate(eff_stress_array, type_facet, gt, facet_nbs);
  Array<Real> dummy_array(
      facet_loop.size() * fe_engine.getNbIntegrationPoints(type_facet), 1, 1.);
  Real loop_area = fe_engine.integrate(dummy_array, type_facet, gt, facet_nbs);
  return int_eff_stress / loop_area;
}

/* --------------------------------------------------------------- */
void RVETools::decomposeFacetLoopPerMaterial(
    const Array<Element> & facet_loop, UInt crack_nb,
    std::map<UInt, std::map<UInt, UInt>> & mat_index_facet_nbs_crack_nbs) {
  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);

  for (auto & facet : facet_loop) {
    ElementType type_facet = facet.type;
    GhostType gt = facet.ghost_type;
    auto & facet_material_index =
        coh_model.getFacetMaterial(type_facet, gt)(facet.element);
    mat_index_facet_nbs_crack_nbs[facet_material_index][facet.element] =
        crack_nb;
  }
}

/* --------------------------------------------------------------- */
void RVETools::applyEigenOpening(Real eigen_strain) {
  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);

  if (not coh_model.getIsExtrinsic()) {
    AKANTU_EXCEPTION(
        "This function can only be used for extrinsic cohesive elements");
  }

  const Mesh & mesh_facets = coh_model.getMeshFacets();
  const UInt dim = model.getSpatialDimension();

  for (auto && mat : model.getMaterials()) {
    auto * mat_coh = dynamic_cast<MaterialCohesive *>(&mat);

    if (mat_coh == nullptr)
      continue;

    for (auto gt : ghost_types) {
      for (auto && type_facet : mesh_facets.elementTypes(dim - 1)) {
        ElementType type_cohesive =
            FEEngine::getCohesiveElementType(type_facet);
        UInt nb_quad_cohesive = model.getFEEngine("CohesiveFEEngine")
                                    .getNbIntegrationPoints(type_cohesive);
        if (not model.getMesh().getNbElement(type_cohesive, gt))
          continue;

        const Array<UInt> & elem_filter =
            mat_coh->getElementFilter(type_cohesive, gt);

        if (not elem_filter.size()) {
          continue;
        }

        // access openings of cohesive elements
        auto & eig_opening = mat_coh->getEigenOpening(type_cohesive, gt);
        auto & normals = mat_coh->getNormals(type_cohesive, gt);
        auto eig_op_it = eig_opening.begin(dim);
        auto normals_it = normals.begin(dim);

        // apply eigen strain as opening
        for (auto el_order : arange(elem_filter.size())) {

          Element cohesive{type_cohesive, elem_filter(el_order), gt};
          const Vector<Element> & facets_to_cohesive =
              mesh_facets.getSubelementToElement(cohesive);
          Real facet_indiam = MeshUtils::getInscribedCircleDiameter(
              model, facets_to_cohesive(0));
          auto eigen_opening = eigen_strain * facet_indiam;

          Vector<Real> normal(normals_it[el_order * nb_quad_cohesive]);
          for (UInt i : arange(nb_quad_cohesive)) {
            auto eig_op = eig_op_it[el_order * nb_quad_cohesive + i];

            eig_op = normal * eigen_opening;
          }
        }
      }
    }
  }
}
/* --------------------------------------------------------------- */
void RVETools::applyEigenOpeningFluidVolumeBased(Real fluid_volume_ratio) {
  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);

  if (not coh_model.getIsExtrinsic()) {
    AKANTU_EXCEPTION(
        "This function can only be used for extrinsic cohesive elements");
  }

  const Mesh & mesh_facets = coh_model.getMeshFacets();
  const UInt dim = model.getSpatialDimension();
  // compute total asr gel volume
  if (!this->volume)
    computeModelVolume();
  Real fluid_volume = this->volume * fluid_volume_ratio;

  // compute total area of existing cracks
  auto data_agg = computeCrackData("agg-agg");
  auto data_agg_mor = computeCrackData("agg-mor");
  auto data_mor = computeCrackData("mor-mor");
  auto area_agg = std::get<0>(data_agg);
  auto area_agg_mor = std::get<0>(data_agg_mor);
  auto area_mor = std::get<0>(data_mor);
  Real total_area = area_agg + area_agg_mor + area_mor;

  // compute necessary opening to cover fit gel volume by crack area
  Real average_opening = fluid_volume / total_area;

  for (auto && mat : model.getMaterials()) {
    auto * mat_coh = dynamic_cast<MaterialCohesive *>(&mat);

    if (mat_coh == nullptr)
      continue;

    for (auto gt : ghost_types) {
      for (auto && type_facet : mesh_facets.elementTypes(dim - 1)) {
        ElementType type_cohesive =
            FEEngine::getCohesiveElementType(type_facet);
        UInt nb_quad_cohesive = model.getFEEngine("CohesiveFEEngine")
                                    .getNbIntegrationPoints(type_cohesive);
        if (not model.getMesh().getNbElement(type_cohesive, gt))
          continue;

        const Array<UInt> & elem_filter =
            mat_coh->getElementFilter(type_cohesive, gt);

        if (not elem_filter.size()) {
          continue;
        }

        // access openings of cohesive elements
        auto & eig_opening = mat_coh->getEigenOpening(type_cohesive, gt);
        auto & normals = mat_coh->getNormals(type_cohesive, gt);
        auto eig_op_it = eig_opening.begin(dim);
        auto normals_it = normals.begin(dim);

        // apply average opening directly as eigen
        for (auto el_order : arange(elem_filter.size())) {
          Vector<Real> normal(normals_it[el_order * nb_quad_cohesive]);
          for (UInt i : arange(nb_quad_cohesive)) {
            auto eig_op = eig_op_it[el_order * nb_quad_cohesive + i];
            eig_op = normal * average_opening;
          }
        }
      }
    }
  }
}
/* --------------------------------------------------------------- */
void RVETools::applyPointForceToCrackCentralNodes(Real force_norm) {
  auto && mesh = model.getMesh();
  auto dim = mesh.getSpatialDimension();
  auto & forces = model.getExternalForce();
  auto & pos = mesh.getNodes();
  auto it_force = make_view(forces, dim).begin();
  auto it_pos = make_view(pos, dim).begin();

  for (auto && data : zip(crack_central_node_pairs, crack_normals_pairs)) {
    auto && node_pair = std::get<0>(data);
    auto && normals_pair = std::get<1>(data);
    auto node1 = node_pair.first;
    auto node2 = node_pair.second;
    auto normal1 = normals_pair.first;
    auto normal2 = normals_pair.second;
    Vector<Real> force1 = it_force[node1];
    Vector<Real> force2 = it_force[node2];
    force1 = normal1 * force_norm / 2.;
    force2 = normal2 * force_norm / 2.;
  }
}
/* --------------------------------------------------------------- */
void RVETools::applyPointForceDelayed(Real loading_rate,
                                      const Array<Real> loading_times,
                                      Real time, Real multiplier) {
  auto && mesh = model.getMesh();
  auto dim = mesh.getSpatialDimension();
  auto & forces = model.getExternalForce();
  auto & pos = mesh.getNodes();
  auto it_force = make_view(forces, dim).begin();
  auto it_pos = make_view(pos, dim).begin();
  AKANTU_DEBUG_ASSERT(loading_times.size() == crack_central_node_pairs.size(),
                      "Number of starting loading times does not correspond to "
                      "the number of node pairs");

  for (auto && data :
       zip(crack_central_node_pairs, crack_normals_pairs, loading_times)) {
    auto && node_pair = std::get<0>(data);
    auto && normals_pair = std::get<1>(data);
    auto && loading_time = std::get<2>(data);
    // skip if the crack should not be loaded YET
    if (time < loading_time)
      continue;
    Real force_norm = loading_rate * (time - loading_time) * multiplier;
    auto node1 = node_pair.first;
    auto node2 = node_pair.second;
    auto normal1 = normals_pair.first;
    auto normal2 = normals_pair.second;
    Vector<Real> force1 = it_force[node1];
    Vector<Real> force2 = it_force[node2];
    force1 = normal1 * force_norm / 2.;
    force2 = normal2 * force_norm / 2.;
  }
}
/* --------------------------------------------------------------- */
void RVETools::applyPointForceDistributed(Real radius, Real F) {
  AKANTU_DEBUG_ASSERT(radius != 0,
                      "Expanding pocket radius could not be equal to 0");
  auto && mesh = model.getMesh();
  auto & mesh_facets = mesh.getMeshFacets();
  auto dim = mesh.getSpatialDimension();
  auto gt = _not_ghost;
  auto & forces = model.getExternalForce();
  auto & pos = mesh.getNodes();
  auto it_force = make_view(forces, dim).begin();
  auto it_pos = make_view(pos, dim).begin();

  for (auto && data : zip(crack_central_node_pairs, crack_normals_pairs)) {
    std::map<std::pair<UInt, UInt>, std::pair<Real, Vector<Real>>>
        loaded_nodes_distance_normal;
    auto central_node_pair = std::get<0>(data);
    auto normals_pair = std::get<1>(data);
    if (central_node_pair.first < central_node_pair.second) {
      loaded_nodes_distance_normal[central_node_pair] =
          std::make_pair(0., normals_pair.first);
    } else {
      auto copy = central_node_pair;
      central_node_pair.first = copy.second;
      central_node_pair.second = copy.first;
      loaded_nodes_distance_normal[central_node_pair] =
          std::make_pair(0., normals_pair.second);
    }
    auto central_node = central_node_pair.first;
    Vector<Real> sphere_center = it_pos[central_node];
    // normal opening of central nodes
    Vector<Real> central_opening = it_pos[central_node_pair.first];
    central_opening -= Vector<Real>(it_pos[central_node_pair.second]);
    auto normal_centr_open_norm = central_opening.dot(normals_pair.first);
    normal_centr_open_norm = std::abs(normal_centr_open_norm);

    // search for cohesives falling into the sphere
    for (auto coh_type : mesh.elementTypes(dim, gt, _ek_cohesive)) {
      auto && coh_conn = mesh.getConnectivity(coh_type, gt);
      auto nb_coh_nodes = coh_conn.getNbComponent();
      auto coh_conn_it = coh_conn.begin(coh_conn.getNbComponent());
      UInt nb_quad_cohesive = model.getFEEngine("CohesiveFEEngine")
                                  .getNbIntegrationPoints(coh_type);

      for (auto && mat : model.getMaterials()) {
        auto * mat_coh = dynamic_cast<MaterialCohesive *>(&mat);

        if (mat_coh == nullptr)
          continue;

        const auto & filter_map = mat.getElementFilter();
        if (!filter_map.exists(coh_type, gt))
          continue;

        const Array<UInt> & filter = mat.getElementFilter(coh_type, gt);
        if (filter.size() == 0)
          continue;

        auto && damage = mat.getInternal<Real>("damage")(coh_type, gt);
        auto && coh_normals = mat_coh->getNormals(coh_type, gt);
        auto coh_norm_it = coh_normals.begin(dim);
        auto && normal_opening_norm =
            mat_coh->getInternal<Real>("normal_opening_norm")(coh_type, gt);

        for (auto && data : enumerate(filter)) {
          auto local_coh_nb = std::get<0>(data);
          auto coh_el_nb = std::get<1>(data);
          // consider only fully broken elements
          if (damage(local_coh_nb * nb_quad_cohesive) < 1.)
            continue;
          // consider only elements with big opening
          if (normal_opening_norm(local_coh_nb * nb_quad_cohesive) <
              0.4 * normal_centr_open_norm)
            continue;
          Element coh_el{coh_type, coh_el_nb, gt};
          auto coh_facet = mesh_facets.getSubelementToElement(coh_el)(0);
          Vector<Real> coh_bary(dim);
          mesh_facets.getBarycenter(coh_facet, coh_bary);
          auto distance = sphere_center.distance(coh_bary);
          // quickly discard if cohesive is further from the center
          if (distance >= radius)
            continue;
          // otherwise check nodes
          auto coh_nodes = coh_conn_it[coh_el_nb];
          std::map<UInt, UInt> counts;

          for (auto i : coh_nodes)
            ++counts[i];
          auto nb_duplicated_nodes =
              std::count_if(counts.begin(), counts.end(),
                            [](auto const & p) { return p.second == 1; });
          // discard cohesives with not all nodes opened
          if (nb_duplicated_nodes != nb_coh_nodes)
            continue;

          // otherwise add all node pairs into the set
          for (UInt i : arange(nb_coh_nodes / 2)) {
            Vector<Real> node_coord = it_pos[coh_nodes(i)];
            auto distance = sphere_center.distance(node_coord);
            // exclude nodes falling outside the sphere
            if (distance <= radius) {
              Vector<Real> coh_normal =
                  coh_norm_it[local_coh_nb * nb_quad_cohesive];
              std::pair<UInt, UInt> nodes_pair;
              Vector<Real> corr_normal(dim);
              if (coh_nodes(i) < coh_nodes(i + nb_coh_nodes / 2)) {
                nodes_pair = std::make_pair(coh_nodes(i),
                                            coh_nodes(i + nb_coh_nodes / 2));
                corr_normal = coh_normal;
              } else {
                nodes_pair = std::make_pair(coh_nodes(i + nb_coh_nodes / 2),
                                            coh_nodes(i));
                corr_normal = -1. * coh_normal;
              }
              if (loaded_nodes_distance_normal.find(nodes_pair) ==
                  loaded_nodes_distance_normal.end()) {
                loaded_nodes_distance_normal[nodes_pair] =
                    std::make_pair(distance, corr_normal);
              }
            }
          }
        }
      }
    }

    // compute central loading
    Real denominator{0};
    for (auto && pair : loaded_nodes_distance_normal) {
      denominator += 1 - pair.second.first / radius;
    }
    auto central_load = F / 2 / denominator;

    // distribute loading between all the node pairs
    for (auto && pair : loaded_nodes_distance_normal) {
      auto node_pair = pair.first;
      auto distance = pair.second.first;
      auto normal = pair.second.second;
      auto node1 = node_pair.first;
      auto node2 = node_pair.second;
      Vector<Real> force1 = it_force[node1];
      Vector<Real> force2 = it_force[node2];
      force1 = normal * central_load * (1 - distance / radius);
      force2 = -1. * normal * central_load * (1 - distance / radius);
    }
  }
}

/* --------------------------------------------------------------- */
void RVETools::applyPointForcesFacetLoop(Real load) {

  auto && mesh = model.getMesh();
  auto dim = mesh.getSpatialDimension();
  auto & forces = model.getExternalForce();
  auto it_force = make_view(forces, dim).begin();

  CSR<Element> nodes_to_cohesives;
  MeshUtils::buildNode2Elements(mesh, nodes_to_cohesives, dim, _ek_cohesive);

  for (auto && data : zip(crack_central_node_pairs)) {
    std::map<std::pair<UInt, UInt>, Vector<Real>> loaded_nodes_normals;
    std::set<std::pair<UInt, UInt>> loaded_nodes;
    auto central_node_pair = std::get<0>(data);
    auto central_node = central_node_pair.first;

    for (auto & coh : nodes_to_cohesives.getRow(central_node)) {
      // check duplicated nodes
      auto coh_nodes = mesh.getConnectivity(coh);
      auto nb_coh_nodes = coh_nodes.size();
      auto nb_quad_cohesive = model.getFEEngine("CohesiveFEEngine")
                                  .getNbIntegrationPoints(coh.type);

      auto mat_id =
          model.getMaterialByElement(coh.type, coh.ghost_type)(coh.element);
      auto coh_loc_num = model.getMaterialLocalNumbering(
          coh.type, coh.ghost_type)(coh.element);
      auto && mat = model.getMaterial(mat_id);
      auto * mat_coh = dynamic_cast<MaterialCohesive *>(&mat);

      auto && coh_normals = mat_coh->getNormals(coh.type, coh.ghost_type);
      auto coh_norm_it = coh_normals.begin(dim);
      Vector<Real> coh_normal = coh_norm_it[coh_loc_num * nb_quad_cohesive];
      std::pair<UInt, UInt> nodes_pair;
      for (UInt i : arange(nb_coh_nodes / 2)) {
        if (coh_nodes(i) != coh_nodes(i + nb_coh_nodes / 2)) {
          nodes_pair.first =
              std::min(coh_nodes(i), coh_nodes(i + nb_coh_nodes / 2));
          nodes_pair.second =
              std::max(coh_nodes(i), coh_nodes(i + nb_coh_nodes / 2));
          Vector<Real> corr_normal = coh_normal;
          if (nodes_pair.first != coh_nodes(i))
            corr_normal *= -1.;
          auto res = loaded_nodes.insert(nodes_pair);
          if (res.second == true) {
            loaded_nodes_normals[nodes_pair] = corr_normal;
          } else {
            loaded_nodes_normals[nodes_pair] += corr_normal;
          }
        }
      }
    }

    // sum up all previously applied forces
    Real current_load{0};
    for (auto && pair : loaded_nodes_normals) {
      auto node_pair = pair.first;
      auto node1 = node_pair.first;
      auto node2 = node_pair.second;
      Vector<Real> force1 = it_force[node1];
      Vector<Real> force2 = it_force[node2];
      current_load += std::abs(force1.norm()) + std::abs(force2.norm());
    }
    Real load_incr = load - current_load;

    // split load increase between all the node pairs
    for (auto && pair : loaded_nodes_normals) {
      auto node_pair = pair.first;
      auto normal = pair.second;
      if (not Math::are_float_equal(normal.norm(), 0))
        normal /= normal.norm();
      auto node1 = node_pair.first;
      auto node2 = node_pair.second;
      Vector<Real> force1 = it_force[node1];
      Vector<Real> force2 = it_force[node2];
      force1 += normal * load_incr / 2. / Real(loaded_nodes_normals.size());
      force2 +=
          -1. * normal * load_incr / 2. / Real(loaded_nodes_normals.size());
    }
  }
}

/* --------------------------------------------------------------- */
void RVETools::outputCrackData(std::ofstream & file_output, Real time) {

  auto data_agg = computeCrackData("agg-agg");
  auto data_agg_mor = computeCrackData("agg-mor");
  auto data_mor = computeCrackData("mor-mor");
  auto area_agg = std::get<0>(data_agg);
  auto vol_agg = std::get<1>(data_agg);
  auto area_agg_mor = std::get<0>(data_agg_mor);
  auto vol_agg_mor = std::get<1>(data_agg_mor);
  auto area_mor = std::get<0>(data_mor);
  auto vol_mor = std::get<1>(data_mor);
  Real total_area = area_agg + area_agg_mor + area_mor;
  Real total_volume = vol_agg + vol_agg_mor + vol_mor;

  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();

  if (prank == 0) {
    file_output << time << "," << area_agg << "," << area_agg_mor << ","
                << area_mor << "," << vol_agg << "," << vol_agg_mor << ","
                << vol_mor << "," << total_area << "," << total_volume
                << std::endl;
  }
}

/* --------------------------------------------------------------- */
std::tuple<Real, Real> RVETools::computeCrackData(const ID & material_name) {
  const auto & mesh = model.getMesh();
  const auto & mesh_facets = mesh.getMeshFacets();
  const auto dim = mesh.getSpatialDimension();
  GhostType gt = _not_ghost;
  Material & mat = model.getMaterial(material_name);
  auto * mat_coh = dynamic_cast<MaterialCohesive *>(&mat);

  AKANTU_DEBUG_ASSERT(mat_coh != nullptr,
                      "Crack data can not be computed on material with id" +
                          material_name);

  const ElementTypeMapArray<UInt> & filter_map = mat_coh->getElementFilter();
  const FEEngine & fe_engine = model.getFEEngine("CohesiveFEEngine");
  Real crack_volume{0};
  Real crack_area{0};

  // Loop over the cohesive element types
  for (auto & element_type : filter_map.elementTypes(dim, gt, _ek_cohesive)) {
    const Array<UInt> & filter = filter_map(element_type);
    if (!filter_map.exists(element_type, gt))
      continue;
    if (filter.size() == 0)
      continue;

    UInt nb_quad_cohesive = model.getFEEngine("CohesiveFEEngine")
                                .getNbIntegrationPoints(element_type);

    // compute normal norm of real opening = opening + eigen opening
    auto opening_it = mat_coh->getOpening(element_type, gt).begin(dim);
    auto eigen_opening = mat_coh->getEigenOpening(element_type, gt);
    auto eigen_opening_it = eigen_opening.begin(dim);
    Array<Real> normal_eigen_opening_norm(
        filter.size() * fe_engine.getNbIntegrationPoints(element_type), 1, 0);
    auto normal_eigen_opening_norm_it = normal_eigen_opening_norm.begin();
    Array<Real> real_normal_opening_norm(
        filter.size() * fe_engine.getNbIntegrationPoints(element_type), 1, 0);
    auto real_normal_opening_norm_it = real_normal_opening_norm.begin();
    auto normal_it = mat_coh->getNormals(element_type, gt).begin(dim);
    auto normal_end = mat_coh->getNormals(element_type, gt).end(dim);
    // loop on each quadrature point
    for (UInt i = 0; normal_it != normal_end; ++i, ++opening_it,
              ++eigen_opening_it, ++normal_it, ++real_normal_opening_norm_it,
              ++normal_eigen_opening_norm_it) {
      UInt coh_nb = filter(floor(i / nb_quad_cohesive));
      Element coh{element_type, coh_nb, gt};
      auto && coh_facets = mesh_facets.getSubelementToElement(coh);
      bool blowed_up{false};
      Real inscr_diam{0};
      for (auto & coh_facet : coh_facets) {
        inscr_diam =
            MeshUtils::getInscribedCircleDiameter(model, coh_facet, true);
        auto & fe_engine = model.getFEEngine("FacetsFEEngine");
        auto initial_area = MeshUtils::getFacetArea(fe_engine, coh_facet);
        auto current_area = MeshUtils::getCurrentFacetArea(model, coh_facet);

        // discard unrealistic deformations (factor is choosen intuitively)
        if (current_area >= 2 * initial_area) {
          blowed_up = true;
          break;
        }
      }

      if (blowed_up)
        continue;

      auto real_opening = *opening_it + *eigen_opening_it;
      auto opening_norm = real_opening.dot(*normal_it);

      // discard unrealistic openings (factor is choosen intuitively)
      if (opening_norm > 10 * inscr_diam)
        continue;
      *real_normal_opening_norm_it = opening_norm;
      Vector<Real> eig_opening(*eigen_opening_it);
      *normal_eigen_opening_norm_it = eig_opening.dot(*normal_it);
    }

    crack_volume =
        fe_engine.integrate(real_normal_opening_norm, element_type, gt, filter);

    Array<Real> area(
        filter.size() * fe_engine.getNbIntegrationPoints(element_type), 1, 1.);
    crack_area = fe_engine.integrate(area, element_type, gt, filter);
  }

  auto && comm = akantu::Communicator::getWorldCommunicator();
  comm.allReduce(crack_volume, SynchronizerOperation::_sum);
  comm.allReduce(crack_area, SynchronizerOperation::_sum);

  return std::make_tuple(crack_area, crack_volume);
}

/* --------------------------------------------------------------- */
void RVETools::outputCrackVolumes(std::ofstream & file_output, Real time) {
  const Communicator & comm = Communicator::getWorldCommunicator();
  UInt prank = comm.whoAmI();
  // shift crack numbers by the number of cracks on lower processors
  UInt tot_nb_cracks{crack_central_node_pairs.size()};
  comm.allReduce(tot_nb_cracks, SynchronizerOperation::_sum);
  Array<Real> crack_volumes(tot_nb_cracks, 1, 0.);

  const FEEngine & fe_engine = model.getFEEngine("CohesiveFEEngine");
  auto && mesh = model.getMesh();
  auto && mesh_facets = mesh.getMeshFacets();
  auto dim = mesh.getSpatialDimension();
  GhostType gt = _not_ghost;

  for (auto coh_type : mesh.elementTypes(dim, gt, _ek_cohesive)) {
    auto & crack_numbers = mesh.getData<UInt>("crack_numbers", coh_type, gt);
    for (auto && mat : model.getMaterials()) {
      auto * mat_coh = dynamic_cast<MaterialCohesive *>(&mat);
      if (mat_coh == nullptr)
        continue;

      const auto & filter_map = mat_coh->getElementFilter();
      if (!filter_map.exists(coh_type, gt))
        continue;

      const Array<UInt> & filter = mat_coh->getElementFilter(coh_type, gt);
      if (filter.size() == 0)
        continue;

      UInt nb_quad_cohesive = model.getFEEngine("CohesiveFEEngine")
                                  .getNbIntegrationPoints(coh_type);

      // compute normal norm of real opening = opening + eigen opening
      auto opening_it = mat_coh->getOpening(coh_type, gt).begin(dim);
      auto eigen_opening = mat_coh->getEigenOpening(coh_type, gt);
      auto eigen_opening_it = eigen_opening.begin(dim);
      Array<Real> normal_eigen_opening_norm(
          filter.size() * fe_engine.getNbIntegrationPoints(coh_type), 1, 0);
      auto normal_eigen_opening_norm_it = normal_eigen_opening_norm.begin();
      Array<Real> real_normal_opening_norm(
          filter.size() * fe_engine.getNbIntegrationPoints(coh_type), 1, 0);
      Array<Real> open_volume(filter.size());
      auto real_normal_opening_norm_it = real_normal_opening_norm.begin();
      auto normal_it = mat_coh->getNormals(coh_type, gt).begin(dim);
      auto normal_end = mat_coh->getNormals(coh_type, gt).end(dim);
      // loop on each quadrature point
      for (UInt i = 0; normal_it != normal_end; ++i, ++opening_it,
                ++eigen_opening_it, ++normal_it, ++real_normal_opening_norm_it,
                ++normal_eigen_opening_norm_it) {
        UInt coh_nb = filter(floor(i / nb_quad_cohesive));
        Element coh{coh_type, coh_nb, gt};
        auto && coh_facets = mesh_facets.getSubelementToElement(coh);
        bool blowed_up{false};
        for (auto & coh_facet : coh_facets) {
          auto inscr_diam_init =
              MeshUtils::getInscribedCircleDiameter(model, coh_facet, true);
          auto inscr_diam_curr =
              MeshUtils::getInscribedCircleDiameter(model, coh_facet, false);
          if (inscr_diam_curr >= 2 * inscr_diam_init) {
            blowed_up = true;
            break;
          }
        }

        if (blowed_up)
          continue;

        auto real_opening = *opening_it + *eigen_opening_it;
        *real_normal_opening_norm_it = real_opening.dot(*normal_it);
        Vector<Real> eig_opening(*eigen_opening_it);
        *normal_eigen_opening_norm_it = eig_opening.dot(*normal_it);
      }

      fe_engine.integrate(real_normal_opening_norm, open_volume, 1, coh_type,
                          gt, filter);
      // distribute individual crack volumes per crack number
      for (UInt i : arange(filter.size())) {
        auto coh_el = filter(i);
        auto crack_nb = crack_numbers(coh_el);
        crack_volumes(crack_nb) += open_volume(i);
      }
    }
  }
  comm.allReduce(crack_volumes, SynchronizerOperation::_sum);
  // output into a file
  if (prank == 0) {
    file_output << time;
    for (auto && single_vol : crack_volumes) {
      file_output << "," << single_vol;
    }
    file_output << std::endl;
  }
}

/* --------------------------------------------------------------- */
void RVETools::outputCrackOpenings(std::ofstream & file_output, Real time) {
  const Communicator & comm = Communicator::getWorldCommunicator();
  UInt prank = comm.whoAmI();
  UInt nb_proc = comm.getNbProc();
  Real global_nb_openings(crack_central_node_pairs.size());
  comm.allReduce(global_nb_openings, SynchronizerOperation::_sum);
  // shift crack numbers by the number of cracks on lower processors
  UInt nb_cracks{crack_central_node_pairs.size()};
  comm.exclusiveScan(nb_cracks);
  UInt max_nb_cracks{crack_central_node_pairs.size()};
  comm.allReduce(max_nb_cracks, SynchronizerOperation::_max);
  Array<Real> openings(nb_proc, max_nb_cracks, -1.);

  auto && mesh = model.getMesh();
  auto dim = mesh.getSpatialDimension();
  auto & disp = model.getDisplacement();
  auto & pos = mesh.getNodes();
  auto it_disp = make_view(disp, dim).begin();
  auto it_pos = make_view(pos, dim).begin();
  auto it_open = make_view(openings, max_nb_cracks).begin();
  auto open_end = make_view(openings, max_nb_cracks).end();

  for (UInt i = 0; i < crack_central_node_pairs.size(); i++) {
    auto && node_pair = crack_central_node_pairs(i);
    auto && normals_pair = crack_normals_pairs(i);

    auto node1 = node_pair.first;
    auto node2 = node_pair.second;
    auto normal1 = normals_pair.first;
    Vector<Real> disp1 = it_disp[node1];
    Vector<Real> disp2 = it_disp[node2];
    Vector<Real> rel_disp = disp1 - disp2;
    Real norm_opening = rel_disp.dot(normal1);
    openings(prank, i) = norm_opening;
  }
  comm.allGather(openings);
  if (prank == 0) {
    file_output << time;
    for (; it_open != open_end; ++it_open) {
      Vector<Real> openings_per_proc = *it_open;
      for (auto && opening : openings_per_proc) {
        if (opening != -1.) {
          file_output << "," << opening;
        }
      }
    }
    file_output << std::endl;
  }
}

/* --------------------------------------------------------------- */
void RVETools::outputCrackOpeningsPerProc(std::ofstream & file_output,
                                          Real time) {

  Array<Real> openings(crack_central_node_pairs.size());
  auto && mesh = model.getMesh();
  auto dim = mesh.getSpatialDimension();
  auto & disp = model.getDisplacement();
  auto & pos = mesh.getNodes();
  auto it_disp = make_view(disp, dim).begin();
  auto it_pos = make_view(pos, dim).begin();

  for (UInt i = 0; i < crack_central_node_pairs.size(); i++) {
    auto && node_pair = crack_central_node_pairs(i);
    auto && normals_pair = crack_normals_pairs(i);

    auto node1 = node_pair.first;
    auto node2 = node_pair.second;
    auto normal1 = normals_pair.first;
    Vector<Real> disp1 = it_disp[node1];
    Vector<Real> disp2 = it_disp[node2];
    Vector<Real> rel_disp = disp1 - disp2;
    Real norm_opening = rel_disp.dot(normal1);
    openings(i) = norm_opening;
  }

  file_output << time;
  for (auto && opening : openings) {
    file_output << "," << opening;
  }
  file_output << std::endl;
}

/* --------------------------------------------------------------- */
void RVETools::identifyCrackCentralNodePairsAndNormals() {
  auto && mesh = model.getMesh();
  auto && mesh_facets = mesh.getMeshFacets();
  auto dim = mesh.getSpatialDimension();
  auto & forces = model.getExternalForce();
  auto & pos = mesh.getNodes();
  auto type_facet = *mesh_facets.elementTypes(dim - 1).begin();
  auto && fe_engine = model.getFEEngine("FacetsFEEngine");
  UInt nb_quad_points = fe_engine.getNbIntegrationPoints(type_facet);

  // get normal to the initial positions
  const auto & init_pos = mesh.getNodes();
  Array<Real> quad_normals(0, dim);
  fe_engine.computeNormalsOnIntegrationPoints(init_pos, quad_normals,
                                              type_facet);
  auto normals_it = quad_normals.begin(dim);
  auto it_force = make_view(forces, dim).begin();
  auto it_pos = make_view(pos, dim).begin();

  this->crack_central_node_pairs.resize(crack_facets.size());
  this->crack_normals_pairs.resize(crack_facets.size());
  for (auto && data :
       zip(crack_facets, crack_central_node_pairs, crack_normals_pairs)) {
    auto && facets_per_site = std::get<0>(data);
    auto & node_pair = std::get<1>(data);
    node_pair = std::make_pair(UInt(-1), UInt(-1));
    auto & crack_normals_pair = std::get<2>(data);

    // find first node
    auto && facet_conn = mesh_facets.getConnectivity(facets_per_site(0));
    for (auto && node : facet_conn) {
      auto ret = crack_central_nodes.find(node);
      if (ret != UInt(-1)) {
        node_pair.first = node;
        break;
      }
    }
    AKANTU_DEBUG_ASSERT(node_pair.first != UInt(-1),
                        "First node was not identified");
    // compute first normal
    Vector<Real> normal(dim, 0);
    for (auto && facet : facets_per_site) {
      // skip out opposite side facets
      auto && facet_nodes = mesh_facets.getConnectivity(facet);
      auto ret =
          std::find(facet_nodes.begin(), facet_nodes.end(), node_pair.first);
      if (ret == facet_nodes.end()) {
        continue;
      }

      // find connected solid element
      auto && facet_neighbors = mesh_facets.getElementToSubelement(facet);
      Element solid_el{ElementNull};
      for (auto && neighbor : facet_neighbors) {
        if (mesh.getKind(neighbor.type) == _ek_regular) {
          solid_el = neighbor;
          break;
        }
      }
      AKANTU_DEBUG_ASSERT(solid_el != ElementNull,
                          "Couldn't identify neighboring solid el");

      // find the out-of-plane solid node
      auto && solid_nodes = mesh.getConnectivity(solid_el);
      auto it =
          std::find(solid_nodes.begin(), solid_nodes.end(), node_pair.first);
      AKANTU_DEBUG_ASSERT(it != solid_nodes.end(),
                          "Solid el does not contain central node");

      UInt out_node(-1);
      for (auto && solid_node : solid_nodes) {
        auto ret =
            std::find(facet_nodes.begin(), facet_nodes.end(), solid_node);
        if (ret == facet_nodes.end()) {
          out_node = solid_node;
          break;
        }
      }
      AKANTU_DEBUG_ASSERT(out_node != UInt(-1),
                          "Couldn't identify out-of-plane node");

      // build the reference vector
      Vector<Real> ref_vector(dim);
      mesh_facets.getBarycenter(facet, ref_vector);
      Vector<Real> pos(mesh.getNodes().begin(dim)[out_node]);
      ref_vector = pos - ref_vector;

      // check if ref vector and the normal are in the same half space
      Vector<Real> facet_normal = normals_it[facet.element * nb_quad_points];
      if (ref_vector.dot(facet_normal) < 0) {
        facet_normal *= -1;
      }
      normal += facet_normal;
    }
    crack_normals_pair.first = normal / normal.norm();
    crack_normals_pair.second = normal / normal.norm() * (-1.);

    // find the second one
    Vector<Real> first_node_pos = it_pos[node_pair.first];
    for (auto && crack_central_node : crack_central_nodes) {
      Vector<Real> node_pos = it_pos[crack_central_node];
      if ((first_node_pos == node_pos) &&
          (crack_central_node != node_pair.first)) {
        node_pair.second = crack_central_node;
        break;
      }
    }
    AKANTU_DEBUG_ASSERT(node_pair.second != UInt(-1),
                        "Second node was not identified");
  }
}

/* ----------------------------------------------------------------- */
std::tuple<std::map<Element, Element>, std::map<Element, UInt>, std::set<UInt>,
           std::map<UInt, std::set<Element>>>
RVETools::determineCrackSurface() {
  AKANTU_DEBUG_IN();

  Mesh & mesh = model.getMesh();
  const Mesh & mesh_facets = mesh.getMeshFacets();
  const UInt dim = mesh.getSpatialDimension();
  std::map<Element, Element> contour_subfacets_coh_el;
  std::map<Element, UInt> surface_subfacets_crack_nb;
  std::set<UInt> surface_nodes;
  std::map<UInt, std::set<Element>> contour_nodes_subfacets;
  std::set<Element> visited_subfacets;

  // loop on cohesives (includes duplicated facets)
  for (auto gt : ghost_types) {
    for (auto & type_cohesive : mesh.elementTypes(dim, gt, _ek_cohesive)) {
      auto nb_cohesives = mesh.getNbElement(type_cohesive, gt);
      if (not nb_cohesives) {
        continue;
      }

      for (auto cohesive_nb : arange(nb_cohesives)) {
        Element cohesive{type_cohesive, cohesive_nb, gt};
        auto && facets_to_cohesive =
            mesh_facets.getSubelementToElement(cohesive);
        for (auto coh_facet : facets_to_cohesive) {
          searchCrackSurfaceInSingleFacet(
              coh_facet, contour_subfacets_coh_el, surface_subfacets_crack_nb,
              surface_nodes, contour_nodes_subfacets, visited_subfacets);
        }
      }
    }
  }

  return std::make_tuple(contour_subfacets_coh_el, surface_subfacets_crack_nb,
                         surface_nodes, contour_nodes_subfacets);
}
/* ----------------------------------------------------------------- */
void RVETools::searchCrackSurfaceInSingleFacet(
    Element & facet, std::map<Element, Element> & contour_subfacets_coh_el,
    std::map<Element, UInt> & /*surface_subfacets_crack_nb*/,
    std::set<UInt> & /*surface_nodes*/,
    std::map<UInt, std::set<Element>> & contour_nodes_subfacets,
    std::set<Element> & visited_subfacets) {

  Mesh & mesh = model.getMesh();
  const Mesh & mesh_facets = mesh.getMeshFacets();

  // iterate through subfacets of a facet searching for cohesives
  auto && subfacets_to_facet = mesh_facets.getSubelementToElement(facet);
  for (auto & subfacet : subfacets_to_facet) {

    // discard subfacets on the border of ghost area
    if (subfacet == ElementNull)
      continue;

    // discard pure ghost subfacets
    bool pure_ghost_subf{false};
    auto && subfacet_conn = mesh_facets.getConnectivity(subfacet);
    for (auto & subfacet_node : subfacet_conn) {
      if (mesh.isPureGhostNode(subfacet_node))
        pure_ghost_subf = true;
    }
    if (pure_ghost_subf)
      continue;

    auto ret = visited_subfacets.emplace(subfacet);
    if (not ret.second)
      continue;

    std::map<Element, UInt> connected_cohesives;
    auto && facets_to_subfacet = mesh_facets.getElementToSubelement(subfacet);
    for (auto & neighbor_facet : facets_to_subfacet) {
      auto & connected_to_facet =
          mesh_facets.getElementToSubelement(neighbor_facet);
      for (auto & connected_el : connected_to_facet) {
        if (connected_el.type == _not_defined)
          goto nextfacet;
        if (mesh.getKind(connected_el.type) == _ek_cohesive) {
          connected_cohesives[connected_el]++;
        }
      }
    nextfacet:;
    }
    if (connected_cohesives.size() == 1 and
        connected_cohesives.begin()->second == 2) {
      contour_subfacets_coh_el[subfacet] = connected_cohesives.begin()->first;
      Vector<UInt> contour_subfacet_nodes =
          mesh_facets.getConnectivity(subfacet);
      for (auto & contour_subfacet_node : contour_subfacet_nodes) {
        contour_nodes_subfacets[contour_subfacet_node].emplace(subfacet);
      }
    }
  }
}

/* ------------------------------------------------------------------- */
template <UInt dim> UInt RVETools::updateDeltaMax() {
  AKANTU_DEBUG_IN();

  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);

  if (not coh_model.getIsExtrinsic()) {
    AKANTU_EXCEPTION(
        "This function can only be used for extrinsic cohesive elements");
  }

  UInt nb_updated_quads(0);
  for (auto && mat : model.getMaterials()) {
    auto * mat_coh =
        dynamic_cast<MaterialCohesiveLinearSequential<dim> *>(&mat);

    if (mat_coh == nullptr)
      continue;
    for (auto type : mat_coh->getElementFilter().elementTypes(dim, _not_ghost,
                                                              _ek_cohesive)) {

      nb_updated_quads += mat_coh->updateDeltaMax(type, _not_ghost);
    }
  }

  // summ all updated quads
  auto && comm = akantu::Communicator::getWorldCommunicator();
  comm.allReduce(nb_updated_quads, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return nb_updated_quads;
}

template UInt RVETools::updateDeltaMax<2>();
template UInt RVETools::updateDeltaMax<3>();

/* ------------------------------------------------------------------- */
template <UInt dim>
std::tuple<Real, Element, UInt> RVETools::getMaxDeltaMaxExcess() {
  AKANTU_DEBUG_IN();

  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  coh_model.interpolateStress();
  coh_model.assembleInternalForces();

  Real proc_max_delta_max_excess{0};
  Element critical_coh_el{ElementNull};
  UInt global_nb_penetr_quads{0};
  for (auto && mat : model.getMaterials()) {
    auto * mat_coh =
        dynamic_cast<MaterialCohesiveLinearSequential<dim> *>(&mat);

    if (mat_coh == nullptr)
      continue;
    for (auto type : mat_coh->getElementFilter().elementTypes(dim, _not_ghost,
                                                              _ek_cohesive)) {
      auto data = mat_coh->computeMaxDeltaMaxExcess(type, _not_ghost);
      auto mat_max_delta_max_excess = std::get<0>(data);
      global_nb_penetr_quads += std::get<2>(data);
      if (mat_max_delta_max_excess > proc_max_delta_max_excess) {
        proc_max_delta_max_excess = mat_max_delta_max_excess;
        critical_coh_el = std::get<1>(data);
      }
    }
  }
  // communicate the highest excess over different processors
  auto global_max_delta_max = proc_max_delta_max_excess;
  auto && comm = akantu::Communicator::getWorldCommunicator();
  comm.allReduce(global_max_delta_max, SynchronizerOperation::_max);

  comm.allReduce(global_nb_penetr_quads, SynchronizerOperation::_sum);

  // find the critical coh between processors
  if (global_max_delta_max) {
    UInt coh_nb = critical_coh_el.element;
    if (proc_max_delta_max_excess != global_max_delta_max) {
      coh_nb = 0;
    }
    comm.allReduce(coh_nb, SynchronizerOperation::_max);
    critical_coh_el.element = coh_nb;
  }

  AKANTU_DEBUG_OUT();
  return std::make_tuple(global_max_delta_max, critical_coh_el,
                         global_nb_penetr_quads);
}

template std::tuple<Real, Element, UInt> RVETools::getMaxDeltaMaxExcess<2>();
template std::tuple<Real, Element, UInt> RVETools::getMaxDeltaMaxExcess<3>();

/* ------------------------------------------------------------------ */
Real RVETools::preventCohesiveInsertionInMaterial(std::string facet_mat_name) {
  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  auto & inserter = coh_model.getElementInserter();
  auto & mesh = model.getMesh();
  auto dim = mesh.getSpatialDimension();
  const GhostType gt = _not_ghost;
  auto facet_material_id = model.getMaterialIndex(facet_mat_name);

  auto & mat = model.getMaterial(facet_material_id);
  auto * mat_coh = dynamic_cast<MaterialCohesive *>(&mat);
  AKANTU_DEBUG_ASSERT(mat_coh != nullptr, "Chosen material is not cohesive");

  auto && facets = mat_coh->getFacetFilter();
  UInt facets_excluded{0};
  for (const auto & el_type : facets.elementTypes(dim - 1, gt, _ek_regular)) {
    auto && element_ids = facets(el_type, gt);
    for (auto && el_id : element_ids) {
      inserter.getCheckFacets(el_type, gt)(el_id) = false;
      facets_excluded++;
    }
  }

  auto && comm = akantu::Communicator::getWorldCommunicator();
  comm.allReduce(facets_excluded, SynchronizerOperation::_sum);
  return facets_excluded;
}

/* ------------------------------------------------------------------- */
template <UInt dim> void RVETools::setUpdateStiffness(bool update_stiffness) {
  AKANTU_DEBUG_IN();

  for (auto && mat : model.getMaterials()) {
    auto * mat_coh =
        dynamic_cast<MaterialCohesiveLinearSequential<dim> *>(&mat);

    if (mat_coh == nullptr)
      continue;

    mat_coh->setUpdateStiffness(update_stiffness);
  }

  AKANTU_DEBUG_OUT();
}

template void RVETools::setUpdateStiffness<2>(bool update_stiffness);
template void RVETools::setUpdateStiffness<3>(bool update_stiffness);

} // namespace akantu
