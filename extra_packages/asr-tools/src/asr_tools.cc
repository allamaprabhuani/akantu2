/**
 * @file   ASR_tools.cc
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @brief  implementation tools for the analysis of ASR samples
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
#include "asr_tools.hh"
#include "aka_voigthelper.hh"
#include "communicator.hh"
#include "material_FE2.hh"
#include "material_iterative_stiffness_reduction.hh"
#include "non_linear_solver.hh"
#include "solid_mechanics_model.hh"
#include <fstream>

/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
ASRTools::ASRTools(SolidMechanicsModel & model)
    : model(model), volume(0.), stress_limit(0), nb_dumps(0),
      doubled_facets_ready(false), doubled_nodes_ready(false), node_pairs(0),
      disp_stored(0, model.getSpatialDimension()),
      ext_force_stored(0, model.getSpatialDimension()),
      boun_stored(0, model.getSpatialDimension()) {

  // register event handler for asr tools
  auto & mesh = model.getMesh();
  // auto const dim = model.getSpatialDimension();
  mesh.registerEventHandler(*this, _ehp_lowest);
  if (mesh.hasMeshFacets())
    mesh.getMeshFacets().registerEventHandler(*this, _ehp_lowest);

  // resize the array of tuples
  node_pairs.resize(2);
  // TODO make it dependent on number of gel sites

  /// find four corner nodes of the RVE
  findCornerNodes();

  // /// resize stress limit according to the dimension size
  // switch (dim) {
  // case 2: {
  //   stress_limit.resize(VoigtHelper<2>::size * VoigtHelper<2>::size);
  //   break;
  // }
  // case 3: {
  //   stress_limit.resize(VoigtHelper<3>::size * VoigtHelper<3>::size);
  //   break;
  // }
  // }
}
/* -------------------------------------------------------------------------- */
void ASRTools::computeModelVolume() {
  auto const dim = model.getSpatialDimension();
  auto & mesh = model.getMesh();
  auto & fem = model.getFEEngine("SolidMechanicsFEEngine");
  GhostType gt = _not_ghost;
  this->volume = 0;
  for (auto element_type : mesh.elementTypes(dim, gt, _ek_not_defined)) {
    Array<Real> Volume(mesh.getNbElement(element_type) *
                           fem.getNbIntegrationPoints(element_type),
                       1, 1.);
    this->volume += fem.integrate(Volume, element_type);
  }
  // auto && comm = akantu::Communicator::getWorldCommunicator();
  // std::cout << "starting reduction";
  // comm.allReduce(this->volume, SynchronizerOperation::_sum);
  // std::cout << " ending reduction" << std::endl;
}
/* ------------------------------------------------------------------------- */
void ASRTools::applyFreeExpansionBC() {
  /// boundary conditions
  const auto & mesh = model.getMesh();
  const Vector<Real> & lowerBounds = mesh.getLowerBounds();
  const Vector<Real> & upperBounds = mesh.getUpperBounds();
  const UInt dim = mesh.getSpatialDimension();
  const Array<Real> & pos = mesh.getNodes();
  Array<Real> & disp = model.getDisplacement();
  Array<bool> & boun = model.getBlockedDOFs();

  /// accessing bounds
  Real bottom = lowerBounds(1);
  Real left = lowerBounds(0);
  Real top = upperBounds(1);
  Real eps = std::abs((top - bottom) * 1e-6);

  switch (dim) {
  case 2: {
    for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
      if (std::abs(pos(i, 1) - bottom) < eps) {
        boun(i, 1) = true;
        disp(i, 1) = 0.0;

        if (std::abs(pos(i, 0) - left) < eps) {
          boun(i, 0) = true;
          disp(i, 0) = 0.0;
        }
      }
    }
    break;
  }
  case 3: {
    /// accessing bounds
    Real back = lowerBounds(2);
    for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
      if (std::abs(pos(i, 1) - bottom) < eps) {
        boun(i, 1) = true;
        disp(i, 1) = 0.0;

        if ((std::abs(pos(i, 0) - left) < eps) &&
            (std::abs(pos(i, 2) - back) < eps)) {
          boun(i, 0) = true;
          boun(i, 2) = true;
          disp(i, 0) = 0.0;
          disp(i, 2) = 0.0;
        }
      }
    }
    break;
  }
  }
}

/* --------------------------------------------------------------------------
 */
void ASRTools::applyLoadedBC(const Vector<Real> & traction,
                             const ID & element_group, bool multi_axial) {

  /// boundary conditions
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
  disp.clear();
  boun.clear();

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

/* --------------------------------------------------------------------------
 */
void ASRTools::fillNodeGroup(NodeGroup & node_group, bool multi_axial) {
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  const Vector<Real> & upperBounds = mesh.getUpperBounds();
  const Vector<Real> & lowerBounds = mesh.getLowerBounds();
  Real top = upperBounds(1);
  Real right = upperBounds(0);
  Real bottom = lowerBounds(1);
  Real front;
  if (dim == 3)
    front = upperBounds(2);
  Real eps = std::abs((top - bottom) * 1e-6);
  const Array<Real> & pos = mesh.getNodes();
  if (multi_axial) {
    /// fill the NodeGroup with the nodes on the left, bottom and back
    /// surface
    for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
      if (std::abs(pos(i, 0) - right) < eps)
        node_group.add(i);
      if (std::abs(pos(i, 1) - top) < eps)
        node_group.add(i);
      if (std::abs(pos(i, 2) - front) < eps)
        node_group.add(i);
    }
  }
  /// fill the NodeGroup with the nodes on the bottom surface
  else {
    for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
      if (std::abs(pos(i, 1) - top) < eps)
        node_group.add(i);
    }
  }
}

/* --------------------------------------------------------------------------
 */
Real ASRTools::computeAverageDisplacement(SpatialDirection direction) {

  Real av_displ = 0;
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  const Vector<Real> & upperBounds = mesh.getUpperBounds();
  const Vector<Real> & lowerBounds = mesh.getLowerBounds();
  Real right = upperBounds(0);
  Real left = lowerBounds(0);
  Real top = upperBounds(1);
  Real front;
  if (dim == 3)
    front = upperBounds(2);
  Real eps = std::abs((right - left) * 1e-6);
  const Array<Real> & pos = mesh.getNodes();
  const Array<Real> & disp = model.getDisplacement();
  UInt nb_nodes = 0;

  if (direction == _x) {
    for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
      if (std::abs(pos(i, 0) - right) < eps && mesh.isLocalOrMasterNode(i)) {
        av_displ += disp(i, 0);
        ++nb_nodes;
      }
    }
  }

  else if (direction == _y) {
    for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
      if (std::abs(pos(i, 1) - top) < eps && mesh.isLocalOrMasterNode(i)) {
        av_displ += disp(i, 1);
        ++nb_nodes;
      }
    }
  }

  else if ((direction == _z) && (model.getSpatialDimension() == 3)) {
    AKANTU_DEBUG_ASSERT(model.getSpatialDimension() == 3,
                        "no z-direction in 2D problem");
    for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
      if (std::abs(pos(i, 2) - front) < eps && mesh.isLocalOrMasterNode(i)) {
        av_displ += disp(i, 2);
        ++nb_nodes;
      }
    }
  }

  else
    AKANTU_EXCEPTION("The parameter for the testing direction is wrong!!!");

  auto && comm = akantu::Communicator::getWorldCommunicator();
  comm.allReduce(av_displ, SynchronizerOperation::_sum);
  comm.allReduce(nb_nodes, SynchronizerOperation::_sum);

  return av_displ / nb_nodes;
}

/* --------------------------------------------------------------------------
 */
Real ASRTools::computeVolumetricExpansion(SpatialDirection direction) {

  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  UInt tensor_component = 0;
  if (dim == 2) {
    if (direction == _x)
      tensor_component = 0;
    else if (direction == _y)
      tensor_component = 3;
    else
      AKANTU_EXCEPTION("The parameter for the testing direction is wrong!!!");
  }

  else if (dim == 3) {
    if (direction == _x)
      tensor_component = 0;
    else if (direction == _y)
      tensor_component = 4;
    else if (direction == _z)
      tensor_component = 8;
    else
      AKANTU_EXCEPTION("The parameter for the testing direction is wrong!!!");
  }

  else
    AKANTU_EXCEPTION("This example does not work for 1D!!!");

  Real gradu_tot = 0;
  Real tot_volume = 0;
  GhostType gt = _not_ghost;

  for (auto element_type : mesh.elementTypes(dim, gt, _ek_not_defined)) {
    const FEEngine & fe_engine = model.getFEEngine();
    for (UInt m = 0; m < model.getNbMaterials(); ++m) {
      const ElementTypeMapUInt & element_filter_map =
          model.getMaterial(m).getElementFilter();
      if (!element_filter_map.exists(element_type, gt))
        continue;
      const Array<UInt> & elem_filter =
          model.getMaterial(m).getElementFilter(element_type);
      if (!elem_filter.size())
        continue;
      // const Array<Real> & gradu_vec =
      // model.getMaterial(m).getGradU(element_type);
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

  auto && comm = akantu::Communicator::getWorldCommunicator();
  comm.allReduce(gradu_tot, SynchronizerOperation::_sum);
  comm.allReduce(tot_volume, SynchronizerOperation::_sum);

  return gradu_tot / tot_volume;
}

/* --------------------------------------------------------------------------
 */
Real ASRTools::computeDamagedVolume(const ID & mat_name) {
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  Real total_damage = 0;
  GhostType gt = _not_ghost;
  const Material & mat = model.getMaterial(mat_name);

  for (auto element_type : mesh.elementTypes(dim, gt, _ek_not_defined)) {
    const FEEngine & fe_engine = model.getFEEngine();
    const ElementTypeMapUInt & element_filter_map = mat.getElementFilter();
    if (!element_filter_map.exists(element_type, gt))
      continue;
    const Array<UInt> & elem_filter = mat.getElementFilter(element_type);
    if (!elem_filter.size())
      continue;
    const Array<Real> & damage = mat.getInternal<Real>("damage")(element_type);

    total_damage += fe_engine.integrate(damage, element_type, gt, elem_filter);
  }

  auto && comm = akantu::Communicator::getWorldCommunicator();
  comm.allReduce(total_damage, SynchronizerOperation::_sum);

  return total_damage;
}

/* --------------------------------------------------------------------------
 */
void ASRTools::computeStiffnessReduction(std::ofstream & file_output, Real time,
                                         bool tension) {

  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  AKANTU_DEBUG_ASSERT(dim != 1, "Example does not work for 1D");

  /// variables for parallel execution
  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();

  /// save nodal values before test
  storeNodalFields();

  if (dim == 2) {
    Real int_residual_x = performLoadingTest(_x, tension);
    Real int_residual_y = performLoadingTest(_y, tension);
    if (prank == 0)
      file_output << time << "," << int_residual_x << "," << int_residual_y
                  << std::endl;
  }

  else {
    Real int_residual_x = performLoadingTest(_x, tension);
    Real int_residual_y = performLoadingTest(_y, tension);
    Real int_residual_z = performLoadingTest(_z, tension);
    if (prank == 0)
      file_output << time << "," << int_residual_x << "," << int_residual_y
                  << "," << int_residual_z << std::endl;
  }

  /// return the nodal values
  restoreNodalFields();
}
/* --------------------------------------------------------------------------
 */
void ASRTools::computeStiffnessReductionFe2Material(std::ofstream & file_output,
                                                    Real time, bool tension) {

  /// variables for parallel execution
  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();

  /// access the FE2 material and enable usage of homogenized stiffness
  const UInt dim = 2;
  MaterialFE2<dim> & mat =
      dynamic_cast<MaterialFE2<dim> &>(model.getMaterial("FE2_mat"));
  mat.enableHomogenStiffness();

  /// save nodal values before test
  storeNodalFields();

  if (dim == 2) {
    Real int_residual_x = performLoadingTest(_x, tension);
    Real int_residual_y = performLoadingTest(_y, tension);
    if (prank == 0)
      file_output << time << "," << int_residual_x << "," << int_residual_y
                  << std::endl;
  }

  else {
    Real int_residual_x = performLoadingTest(_x, tension);
    Real int_residual_y = performLoadingTest(_y, tension);
    Real int_residual_z = performLoadingTest(_z, tension);
    if (prank == 0)
      file_output << time << "," << int_residual_x << "," << int_residual_y
                  << "," << int_residual_z << std::endl;
  }

  /// return the nodal values
  restoreNodalFields();

  /// disable usage of homogenized stiffness
  mat.disableHomogenStiffness();
}

/* -------------------------------------------------------------------------- */
void ASRTools::storeNodalFields() {
  auto & disp = this->model.getDisplacement();
  auto & boun = this->model.getBlockedDOFs();
  auto & ext_force = this->model.getExternalForce();
  this->disp_stored.copy(disp);
  this->boun_stored.copy(boun);
  this->ext_force_stored.copy(ext_force);
}
/* -------------------------------------------------------------------------- */
void ASRTools::restoreNodalFields() {
  auto & disp = this->model.getDisplacement();
  auto & boun = this->model.getBlockedDOFs();
  auto & ext_force = this->model.getExternalForce();
  disp.copy(this->disp_stored);
  boun.copy(this->boun_stored);
  ext_force.copy(this->ext_force_stored);
  /// update grad_u
  // model.assembleInternalForces();
}

/* ---------------------------------------------------------------------- */
void ASRTools::resetNodalFields() {
  auto & disp = this->model.getDisplacement();
  auto & boun = this->model.getBlockedDOFs();
  auto & ext_force = this->model.getExternalForce();
  disp.clear();
  boun.clear();
  ext_force.clear();
}

/* -------------------------------------------------------------------------- */
void ASRTools::restoreInternalFields() {
  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    Material & mat = model.getMaterial(m);
    mat.restorePreviousState();
  }
}
/* -------------------------------------------------------------------------- */
void ASRTools::resetInternalFields() {
  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    Material & mat = model.getMaterial(m);
    mat.resetInternalsWithHistory();
  }
}
/* --------------------------------------------------------------------------
 */
Real ASRTools::performLoadingTest(SpatialDirection direction, bool tension) {
  UInt dir;

  if (direction == _x)
    dir = 0;
  if (direction == _y)
    dir = 1;
  if (direction == _z) {
    AKANTU_DEBUG_ASSERT(model.getSpatialDimension() == 3,
                        "Error in problem dimension!!!");
    dir = 2;
  }

  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();

  /// boundary conditions
  const Vector<Real> & lowerBounds = mesh.getLowerBounds();
  const Vector<Real> & upperBounds = mesh.getUpperBounds();
  Real bottom = lowerBounds(1);
  Real top = upperBounds(1);
  Real left = lowerBounds(0);
  Real back;

  if (dim == 3)
    back = lowerBounds(2);

  Real eps = std::abs((top - bottom) * 1e-6);
  const auto & pos = mesh.getNodes();
  auto & disp = model.getDisplacement();
  auto & boun = model.getBlockedDOFs();
  auto & ext_force = model.getExternalForce();
  UInt nb_nodes = mesh.getNbNodes();

  disp.clear();
  boun.clear();
  ext_force.clear();

  if (dim == 3) {
    for (UInt i = 0; i < nb_nodes; ++i) {

      /// fix one corner node to avoid sliding
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
        disp(i, dir) = (2 * tension - 1) * 1.e-1;
      }
    }
  } else {
    for (UInt i = 0; i < nb_nodes; ++i) {

      /// fix one corner node to avoid sliding
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
        disp(i, dir) = (2 * tension - 1) * 1.e-1;
      }
    }
  }

  try {
    model.solveStep();
  } catch (debug::Exception & e) {
    auto & solver = model.getNonLinearSolver("static");
    int nb_iter = solver.get("nb_iterations");
    std::cout << "Didn't converge in " << nb_iter << " iterations."
              << std::endl;
    throw e;
  }

  /// compute the force (residual in this case) along the edge of the
  /// imposed displacement
  Real int_residual = 0.;
  const Array<Real> & residual = model.getInternalForce();

  for (UInt n = 0; n < nb_nodes; ++n) {
    if (std::abs(pos(n, dir) - upperBounds(dir)) < eps &&
        mesh.isLocalOrMasterNode(n)) {
      int_residual += -residual(n, dir);
    }
  }

  /// restore historical internal fields
  restoreInternalFields();

  auto && comm = akantu::Communicator::getWorldCommunicator();
  comm.allReduce(int_residual, SynchronizerOperation::_sum);

  return int_residual;
}

/* --------------------------------------------------------------------------
 */
void ASRTools::computeAverageProperties(std::ofstream & file_output) {

  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  AKANTU_DEBUG_ASSERT(dim != 1, "Example does not work for 1D");

  /// variables for parallel execution
  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();

  Real av_strain_x = computeVolumetricExpansion(_x);
  Real av_strain_y = computeVolumetricExpansion(_y);
  Real av_displ_x = computeAverageDisplacement(_x);
  Real av_displ_y = computeAverageDisplacement(_y);
  Real damage_agg = computeDamagedVolume("aggregate");
  Real damage_paste = computeDamagedVolume("paste");

  if (dim == 2) {

    if (prank == 0)
      file_output << av_strain_x << " " << av_strain_y << " " << av_displ_x
                  << " " << av_displ_y << " " << damage_agg << " "
                  << damage_paste << std::endl;
  }

  else {
    Real av_displ_z = computeAverageDisplacement(_z);
    Real av_strain_z = computeVolumetricExpansion(_z);

    if (prank == 0)
      file_output << av_strain_x << " " << av_strain_y << " " << av_strain_z
                  << " " << av_displ_x << " " << av_displ_y << " " << av_displ_z
                  << " " << damage_agg << " " << damage_paste << std::endl;
  }
}

/* --------------------------------------------------------------------------
 */
void ASRTools::computeAverageProperties(std::ofstream & file_output,
                                        Real time) {
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();

  AKANTU_DEBUG_ASSERT(dim != 1, "Example does not work for 1D");

  /// variables for parallel execution
  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();

  Real av_strain_x = computeVolumetricExpansion(_x);
  Real av_strain_y = computeVolumetricExpansion(_y);
  Real av_displ_x = computeAverageDisplacement(_x);
  Real av_displ_y = computeAverageDisplacement(_y);
  Real damage_agg = computeDamagedVolume("aggregate");
  Real damage_paste = computeDamagedVolume("paste");

  if (dim == 2) {

    if (prank == 0)
      file_output << time << "," << av_strain_x << "," << av_strain_y << ","
                  << av_displ_x << "," << av_displ_y << "," << damage_agg << ","
                  << damage_paste << std::endl;
  }

  else {
    Real av_displ_z = computeAverageDisplacement(_z);
    Real av_strain_z = computeVolumetricExpansion(_z);

    if (prank == 0)
      file_output << time << "," << av_strain_x << "," << av_strain_y << ","
                  << av_strain_z << "," << av_displ_x << "," << av_displ_y
                  << "," << av_displ_z << "," << damage_agg << ","
                  << damage_paste << std::endl;
  }
}
/* --------------------------------------------------------------------------
 */
void ASRTools::computeAveragePropertiesFe2Material(std::ofstream & file_output,
                                                   Real time) {
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();

  AKANTU_DEBUG_ASSERT(dim != 1, "Example does not work for 1D");

  /// variables for parallel execution
  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();

  Real av_strain_x = computeVolumetricExpansion(_x);
  Real av_strain_y = computeVolumetricExpansion(_y);
  Real av_displ_x = computeAverageDisplacement(_x);
  Real av_displ_y = computeAverageDisplacement(_y);
  Real damage_paste = averageScalarField("damage_ratio_paste");
  Real damage_agg = averageScalarField("damage_ratio_agg");

  if (dim == 2) {

    if (prank == 0)
      file_output << time << "," << av_strain_x << "," << av_strain_y << ","
                  << av_displ_x << "," << av_displ_y << "," << damage_agg << ","
                  << damage_paste << std::endl;
  }

  else {
    Real av_displ_z = computeAverageDisplacement(_z);
    Real av_strain_z = computeVolumetricExpansion(_z);

    if (prank == 0)
      file_output << time << "," << av_strain_x << "," << av_strain_y << ","
                  << av_strain_z << "," << av_displ_x << "," << av_displ_y
                  << "," << av_displ_z << "," << damage_agg << ","
                  << damage_paste << std::endl;
  }
}

/* --------------------------------------------------------------------------
 */
void ASRTools::saveState(UInt current_step) {

  /// variables for parallel execution
  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();
  GhostType gt = _not_ghost;

  /// save the reduction step on each quad
  UInt nb_materials = model.getNbMaterials();

  /// open the output file
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
    /// get the reduction steps internal field
    const InternalField<UInt> & reduction_steps =
        mat.getInternal<UInt>("damage_step");
    /// loop over all types in that material
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

  /// close the file
  file_output.close();
  comm.barrier();
  /// write the number of the last successfully saved step in a file
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

/* --------------------------------------------------------------------------
 */
bool ASRTools::loadState(UInt & current_step) {
  current_step = 0;
  bool restart = false;

  /// variables for parallel execution
  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();
  GhostType gt = _not_ghost;

  /// proc 0 has to save the current step
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

  /// save the reduction step on each quad
  UInt nb_materials = model.getNbMaterials();

  /// open the output file
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
    /// get the material parameters
    Real E = mat.getParam("E");
    Real max_damage = mat.getParam("max_damage");
    Real max_reductions = mat.getParam("max_reductions");

    /// get the internal field that need to be set
    InternalField<UInt> & reduction_steps =
        const_cast<InternalField<UInt> &>(mat.getInternal<UInt>("damage_step"));
    InternalField<Real> & Sc =
        const_cast<InternalField<Real> &>(mat.getInternal<Real>("Sc"));

    InternalField<Real> & damage =
        const_cast<InternalField<Real> &>(mat.getInternal<Real>("damage"));
    /// loop over all types in that material
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
  /// close the file
  file_input.close();
  return restart;
}

/* --------------------------------------------------------------------------
 */
Real ASRTools::computePhaseVolume(const ID & mat_name) {
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

  auto && comm = akantu::Communicator::getWorldCommunicator();
  comm.allReduce(total_volume, SynchronizerOperation::_sum);

  return total_volume;
}

/* --------------------------------------------------------------------------
 */
void ASRTools::applyEigenGradUinCracks(
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
void ASRTools::applyEigenGradUinCracks(
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

// /*
// --------------------------------------------------------------------------
// */ void ASRTools<dim>::fillCracks(ElementTypeMapReal & saved_damage) {
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

// /*
// --------------------------------------------------------------------------
// */ void ASRTools<dim>::drainCracks(const ElementTypeMapReal &
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
template <UInt dim> Real ASRTools::computeSmallestElementSize() {

  const auto & mesh = model.getMesh();
  //  constexpr auto dim = model.getSpatialDimension();
  /// compute smallest element size
  const Array<Real> & pos = mesh.getNodes();
  Real el_h_min = 10000.;
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

  /// find global Gauss point with highest stress
  auto && comm = akantu::Communicator::getWorldCommunicator();
  comm.allReduce(el_h_min, SynchronizerOperation::_min);

  return el_h_min;
}

// /*
// --------------------------------------------------------------------------
// */
// /// apply homogeneous temperature on the whole specimen

// void ASRTools<dim>::applyTemperatureFieldToSolidmechanicsModel(
//     const Real & temperature) {

//   UInt dim = model.getMesh().getSpatialDimension();

//   for (UInt m = 0; m < model.getNbMaterials(); ++m) {
//     Material & mat = model.getMaterial(m);

//     for (auto el_type : mat.getElementFilter().elementTypes(dim)) {

//       const auto & filter = mat.getElementFilter()(el_type);
//       if (filter.size() == 0)
//         continue;

//       auto & delta_T = mat.getArray<Real>("delta_T", el_type);
//       auto dt = delta_T.begin();
//       auto dt_end = delta_T.end();

//       for (; dt != dt_end; ++dt) {
//         *dt = temperature;
//       }
//     }
//   }
// }

/* --------------------------------------------------------------------------
 */
Real ASRTools::computeDeltaGelStrainThermal(const Real delta_time_day,
                                            const Real k,
                                            const Real activ_energy,
                                            const Real R, const Real T,
                                            Real & amount_reactive_particles,
                                            const Real saturation_const) {
  /// compute increase in gel strain value for interval of time delta_time
  /// as temperatures are stored in C, conversion to K is done

  Real delta_strain = amount_reactive_particles * k *
                      std::exp(-activ_energy / (R * (T + 273.15))) *
                      delta_time_day;

  amount_reactive_particles -= std::exp(-activ_energy / (R * (T + 273.15))) *
                               delta_time_day / saturation_const;

  if (amount_reactive_particles < 0.)
    amount_reactive_particles = 0.;
  return delta_strain;
}

/* -----------------------------------------------------------------------*/
Real ASRTools::computeDeltaGelStrainLinear(const Real delta_time,
                                           const Real k) {
  /// compute increase in gel strain value for dt simply by deps = k *
  /// delta_time

  Real delta_strain = k * delta_time;

  return delta_strain;
}

/* ---------------------------------------------------------------------- */
void ASRTools::applyBoundaryConditionsRve(
    const Matrix<Real> & displacement_gradient) {
  AKANTU_DEBUG_IN();

  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  // /// transform into Green strain to exclude rotations
  // auto F = displacement_gradient + Matrix<Real>::eye(dim);
  // auto E = 0.5 * (F.transpose() * F - Matrix<Real>::eye(dim));
  /// get the position of the nodes
  const Array<Real> & pos = mesh.getNodes();
  /// storage for the coordinates of a given node and the displacement
  /// that will be applied
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
void ASRTools::applyHomogeneousTemperature(const Real & temperature) {
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
void ASRTools::removeTemperature() {
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
void ASRTools::findCornerNodes() {
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
      // node 1
      if (Math::are_float_equal(distance, 0)) {
        corner_nodes(0) = node;
      }
      // node 2
      else if (Math::are_float_equal(X(_x), upper_bounds(_x)) &&
               Math::are_float_equal(X(_y), lower_bounds(_y))) {
        corner_nodes(1) = node;
      }
      // node 3
      else if (Math::are_float_equal(X(_x), upper_bounds(_x)) &&
               Math::are_float_equal(X(_y), upper_bounds(_y))) {
        corner_nodes(2) = node;
      }
      // node 4
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

  //  AKANTU_DEBUG_ASSERT(dim == 2, "This is 2D only!");

  // for (UInt i = 0; i < corner_nodes.size(); ++i) {
  //   if (corner_nodes(i) == UInt(-1))
  //     AKANTU_ERROR("The corner node " << i + 1 << " wasn't found");
  // }
  AKANTU_DEBUG_OUT();
}
/* --------------------------------------------------------------------------
 */
Real ASRTools::averageScalarField(const ID & field_name) {
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

  /// compute total model volume
  if (!this->volume)
    computeModelVolume();

  return average / this->volume;

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
Real ASRTools::averageTensorField(UInt row_index, UInt col_index,
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
          average +=
              int_stress_vec(k, row_index * dim + col_index); // 3 is the value
                                                              // for the yy (in
                                                              // 3D, the value
                                                              // is 4)
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
          /// averaging is done only for normal components, so stress and
          /// strain are equal
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
          /// averaging is done only for normal components, so stress and
          /// strain are equal
          average += int_eigen_gradu_vec(k, row_index * dim + col_index);
      }
    } else {
      AKANTU_ERROR("Averaging not implemented for this field!!!");
    }
  }

  /// compute total model volume
  if (!this->volume)
    computeModelVolume();

  return average / this->volume;

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
void ASRTools::homogenizeStressField(Matrix<Real> & stress) {
  AKANTU_DEBUG_IN();
  stress(0, 0) = averageTensorField(0, 0, "stress");
  stress(1, 1) = averageTensorField(1, 1, "stress");
  stress(0, 1) = averageTensorField(0, 1, "stress");
  stress(1, 0) = averageTensorField(1, 0, "stress");
  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
bool ASRTools::isStressStateTensile(Matrix<Real> & stress) {
  AKANTU_DEBUG_IN();
  Vector<Real> eigenvalues(2);
  stress.eig(eigenvalues);
  Real hydrostatic_stress = 0.5 * (eigenvalues(0) + eigenvalues(1));
  if (hydrostatic_stress > 0)
    return true;
  return false;

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
void ASRTools::homogenizeStiffness(Matrix<Real> & C_macro, bool tensile_test) {
  AKANTU_DEBUG_IN();
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();

  AKANTU_DEBUG_ASSERT(dim == 2, "Is only implemented for 2D!!!");

  /// apply three independent loading states to determine C
  /// 1. eps_el = (1;0;0) 2. eps_el = (0,1,0) 3. eps_el = (0,0,0.5)

  /// clear the eigenstrain TODO why????????
  // Matrix<Real> zero_eigengradu(dim, dim, 0.);
  // GhostType gt = _not_ghost;
  // for (auto element_type : mesh.elementTypes(dim, gt, _ek_not_defined)) {
  //   auto & prestrain_vect =
  //       const_cast<Array<Real>
  //       &>(model.getMaterial("gel").getInternal<Real>(
  //           "eigen_grad_u")(element_type));
  //   auto prestrain_it = prestrain_vect.begin(dim, dim);
  //   auto prestrain_end = prestrain_vect.end(dim, dim);

  //   for (; prestrain_it != prestrain_end; ++prestrain_it)
  //     (*prestrain_it) = zero_eigengradu;
  // }

  /// save nodal values before tests
  storeNodalFields();

  /// storage for results of 3 different loading states
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

  /// save the damage state before filling up cracks
  // ElementTypeMapReal saved_damage("saved_damage");
  // saved_damage.initialize(getFEEngine(), _nb_component = 1,
  // _default_value = 0); this->fillCracks(saved_damage);

  /// virtual test 1:
  H(0, 0) = 0.01 * (2 * tensile_test - 1);
  performVirtualTesting(H, stresses, strains, 0);

  /// virtual test 2:
  H.clear();
  H(1, 1) = 0.01 * (2 * tensile_test - 1);
  performVirtualTesting(H, stresses, strains, 1);

  /// virtual test 3:
  H.clear();
  H(0, 1) = 0.01;
  H(1, 0) = 0.01;
  performVirtualTesting(H, stresses, strains, 2);

  // /// set up the stress limit at 10% of stresses in undamaged state
  // if (first_time) {
  //   auto && str_lim_mat_it =
  //       make_view(this->stress_limit, voigt_size, voigt_size).begin();
  //   *str_lim_mat_it = stresses * 0.1;
  // }

  // /// compare stresses with the lower limit and update if needed
  // auto && str_lim_vec_it = make_view(this->stress_limit,
  // voigt_size).begin(); for (UInt i = 0; i != voigt_size; ++i,
  // ++str_lim_vec_it) {
  //   Vector<Real> stress = stresses(i);
  //   // Vector<Real> stress_lim = stress_limit(i);
  //   Real stress_norm = stress.norm();
  //   Real stress_limit_norm = (*str_lim_vec_it).norm();
  //   if (stress_norm < stress_limit_norm)
  //     // stresses(i) = *str_lim_vec_it;
  //     return;
  // }
  /// drain cracks
  // this->drainCracks(saved_damage);

  /// compute effective stiffness
  Matrix<Real> eps_inverse(voigt_size, voigt_size);
  eps_inverse.inverse(strains);

  /// Make C matrix symmetric
  Matrix<Real> C_direct(voigt_size, voigt_size);
  C_direct.mul<false, false>(stresses, eps_inverse);
  for (UInt i = 0; i != voigt_size; ++i) {
    for (UInt j = 0; j != voigt_size; ++j) {
      C_macro(i, j) = 0.5 * (C_direct(i, j) + C_direct(j, i));
      C_macro(j, i) = C_macro(i, j);
    }
  }

  /// return the nodal values
  restoreNodalFields();

  AKANTU_DEBUG_OUT();
} // namespace akantu

/* --------------------------------------------------------------------------
 */
void ASRTools::performVirtualTesting(const Matrix<Real> & H,
                                     Matrix<Real> & eff_stresses,
                                     Matrix<Real> & eff_strains,
                                     const UInt test_no) {
  AKANTU_DEBUG_IN();
  auto & disp = model.getDisplacement();
  auto & boun = model.getBlockedDOFs();
  auto & ext_force = model.getExternalForce();
  disp.clear();
  boun.clear();
  ext_force.clear();

  applyBoundaryConditionsRve(H);

  auto & solver = model.getNonLinearSolver();
  solver.set("max_iterations", 50);
  solver.set("threshold", 1e-6);
  solver.set("convergence_type", SolveConvergenceCriteria::_solution);

  model.solveStep();

  /// get average stress and strain
  eff_stresses(0, test_no) = averageTensorField(0, 0, "stress");
  eff_strains(0, test_no) = averageTensorField(0, 0, "strain");
  eff_stresses(1, test_no) = averageTensorField(1, 1, "stress");
  eff_strains(1, test_no) = averageTensorField(1, 1, "strain");
  eff_stresses(2, test_no) = averageTensorField(1, 0, "stress");
  eff_strains(2, test_no) = 2. * averageTensorField(1, 0, "strain");

  /// restore historical internal fields
  restoreInternalFields();

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------- */
void ASRTools::homogenizeEigenGradU(Matrix<Real> & eigen_gradu_macro) {
  AKANTU_DEBUG_IN();
  eigen_gradu_macro(0, 0) = averageTensorField(0, 0, "eigen_grad_u");
  eigen_gradu_macro(1, 1) = averageTensorField(1, 1, "eigen_grad_u");
  eigen_gradu_macro(0, 1) = averageTensorField(0, 1, "eigen_grad_u");
  eigen_gradu_macro(1, 0) = averageTensorField(1, 0, "eigen_grad_u");
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------
 */
void ASRTools::fillCracks(ElementTypeMapReal & saved_damage) {
  const auto & mat_gel = model.getMaterial("gel");
  Real E_gel = mat_gel.get("E");
  Real E_homogenized = 0.;

  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    Material & mat = model.getMaterial(m);
    if (mat.getName() == "gel" || mat.getName() == "FE2_mat")
      continue;

    Real E = mat.get("E");
    auto & damage = mat.getInternal<Real>("damage");

    GhostType gt = _not_ghost;
    for (auto && type : damage.elementTypes(gt)) {
      const auto & elem_filter = mat.getElementFilter(type);
      auto nb_integration_point =
          model.getFEEngine().getNbIntegrationPoints(type);

      auto sav_dam_it =
          make_view(saved_damage(type), nb_integration_point).begin();
      for (auto && data :
           zip(elem_filter, make_view(damage(type), nb_integration_point))) {
        auto el = std::get<0>(data);
        auto & dam = std::get<1>(data);
        Vector<Real> sav_dam = sav_dam_it[el];

        sav_dam = dam;

        for (auto q : arange(dam.size())) {
          E_homogenized = (E_gel - E) * dam(q) + E;
          dam(q) = 1. - (E_homogenized / E);
        }
      }
    }
  }
}

/* --------------------------------------------------------------------------
 */
void ASRTools::drainCracks(const ElementTypeMapReal & saved_damage) {
  for (UInt m = 0; m < model.getNbMaterials(); ++m) {
    Material & mat = model.getMaterial(m);
    if (mat.getName() == "gel" || mat.getName() == "FE2_mat")
      continue;
    auto & damage = mat.getInternal<Real>("damage");

    GhostType gt = _not_ghost;
    for (auto && type : damage.elementTypes(gt)) {
      const auto & elem_filter = mat.getElementFilter(type);
      auto nb_integration_point =
          model.getFEEngine().getNbIntegrationPoints(type);

      auto sav_dam_it =
          make_view(saved_damage(type), nb_integration_point).begin();
      for (auto && data :
           zip(elem_filter, make_view(damage(type), nb_integration_point))) {
        auto el = std::get<0>(data);
        auto & dam = std::get<1>(data);
        Vector<Real> sav_dam = sav_dam_it[el];

        dam = sav_dam;
      }
    }
  }
}

/* --------------------------------------------------------------------------
 */
void ASRTools::computeDamageRatio(Real & damage_ratio) {
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

  /// compute total model volume
  if (!this->volume)
    computeModelVolume();

  damage_ratio /= this->volume;
}

/* --------------------------------------------------------------------------
 */
void ASRTools::computeDamageRatioPerMaterial(Real & damage_ratio,
                                             const ID & material_name) {
  const auto & mesh = model.getMesh();
  const auto dim = mesh.getSpatialDimension();
  GhostType gt = _not_ghost;
  Material & mat = model.getMaterial(material_name);
  const ElementTypeMapArray<UInt> & filter_map = mat.getElementFilter();
  const FEEngine & fe_engine = model.getFEEngine();
  damage_ratio = 0.;
  Real mat_volume = 0.;

  // Loop over the boundary element types
  for (auto & element_type : filter_map.elementTypes(dim, gt)) {
    const Array<UInt> & filter = filter_map(element_type);
    if (!filter_map.exists(element_type, gt))
      continue;
    if (filter.size() == 0)
      continue;

    auto & damage_array = mat.getInternal<Real>("damage")(element_type);

    damage_ratio += fe_engine.integrate(damage_array, element_type, gt, filter);
    Array<Real> volume(
        filter.size() * fe_engine.getNbIntegrationPoints(element_type), 1, 1.);
    mat_volume += fe_engine.integrate(volume, element_type, gt, filter);
  }

  damage_ratio /= mat_volume;
}

/* -------------------------------------------------------------------------*/
void ASRTools::dumpRve() {
  //  if (this->nb_dumps % 10 == 0) {
  model.dump();
  //  }
  // this->nb_dumps += 1;
}

/* --------------------------------------------------------------------------
 */
void ASRTools::applyBodyForce() {
  auto spatial_dimension = model.getSpatialDimension();
  model.assembleMassLumped();
  auto & mass = model.getMass();
  auto & force = model.getExternalForce();
  Vector<Real> gravity(spatial_dimension);
  gravity(1) = -9.81;

  for (auto && data : zip(make_view(mass, spatial_dimension),
                          make_view(force, spatial_dimension))) {

    const auto & mass_vec = (std::get<0>(data));
    auto & force_vec = (std::get<1>(data));

    force_vec += gravity * mass_vec;
  }
}
/* -------------------------------------------------------------------------
 */
void ASRTools::insertCohElemByCoord(const Vector<Real> & position) {
  AKANTU_DEBUG_IN();

  auto & mesh = model.getMesh();
  auto dim = mesh.getSpatialDimension();
  const auto gt = _not_ghost;

  // insert cohesive elements
  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  auto & inserter = coh_model.getElementInserter();
  auto & insertion = inserter.getInsertionFacetsByElement();
  auto & mesh_facets = inserter.getMeshFacets();
  Vector<Real> bary_facet(dim);

  Real min_dist = std::numeric_limits<Real>::max();
  Element source_facet;
  for_each_element(mesh_facets,
                   [&](auto && facet) {
                     mesh_facets.getBarycenter(facet, bary_facet);
                     auto dist = bary_facet.distance(position);
                     if (dist < min_dist) {
                       min_dist = dist;
                       source_facet = facet;
                     }
                   },
                   _spatial_dimension = dim - 1);
  inserter.getCheckFacets(source_facet.type, gt)(source_facet.element) = false;
  insertion(source_facet) = true;

  inserter.insertElements();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------
 */
void ASRTools::insertCohElemByLimits(const Matrix<Real> & insertion_limits,
                                     std::string coh_mat_name) {
  AKANTU_DEBUG_IN();
  auto & mesh = model.getMesh();
  auto dim = mesh.getSpatialDimension();
  auto pos_it = make_view(mesh.getNodes(), dim).begin();

  // insert cohesive elements
  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  auto & inserter = coh_model.getElementInserter();
  auto & insertion = inserter.getInsertionFacetsByElement();
  auto & mesh_facets = inserter.getMeshFacets();
  auto coh_material = model.getMaterialIndex(coh_mat_name);
  auto tolerance = Math::getTolerance();
  Vector<Real> bary_facet(dim);

  for_each_element(
      mesh_facets,
      [&](auto && facet) {
        mesh_facets.getBarycenter(facet, bary_facet);
        UInt coord_in_limit = 0;

        while (coord_in_limit < dim and
               bary_facet(coord_in_limit) >
                   (insertion_limits(coord_in_limit, 0) - tolerance) and
               bary_facet(coord_in_limit) <
                   (insertion_limits(coord_in_limit, 1) + tolerance))
          ++coord_in_limit;

        if (coord_in_limit == dim) {
          coh_model.getFacetMaterial(facet.type, facet.ghost_type)(
              facet.element) = coh_material;
          inserter.getCheckFacets(facet.type, facet.ghost_type)(facet.element) =
              false;
          insertion(facet) = true;
        }
      },
      _spatial_dimension = dim - 1);

  inserter.insertElements();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------
 */
void ASRTools::insertCohElemRandomly(const UInt & nb_coh_elem,
                                     std::string coh_mat_name,
                                     std::string matrix_mat_name) {
  AKANTU_DEBUG_IN();

  auto & mesh = model.getMesh();
  auto dim = mesh.getSpatialDimension();
  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);

  // get element types
  const GhostType gt = _not_ghost;
  auto type = *mesh.elementTypes(dim, gt, _ek_regular).begin();
  auto type_facets = Mesh::getFacetType(type);
  auto & inserter = coh_model.getElementInserter();
  auto & insertion = inserter.getInsertionFacetsByElement();
  Vector<Real> bary_facet(dim);
  auto coh_material = model.getMaterialIndex(coh_mat_name);
  auto & mesh_facets = inserter.getMeshFacets();
  auto matrix_material_id = model.getMaterialIndex(matrix_mat_name);

  auto & matrix_elements =
      mesh_facets.createElementGroup("matrix_facets", dim - 1);

  for_each_element(mesh_facets,
                   [&](auto && facet) {
                     mesh_facets.getBarycenter(facet, bary_facet);
                     auto & facet_material = coh_model.getFacetMaterial(
                         facet.type, facet.ghost_type)(facet.element);

                     if (facet_material == matrix_material_id)
                       matrix_elements.add(facet);
                   },
                   _spatial_dimension = dim - 1);

  // Will be used to obtain a seed for the random number engine
  std::random_device rd;
  // Standard mersenne_twister_engine seeded with rd()
  std::mt19937 random_generator(rd());
  std::uniform_int_distribution<> dis(
      0, matrix_elements.getElements(type_facets).size() - 1);

  for (auto _[[gnu::unused]] : arange(nb_coh_elem)) {
    auto id = dis(random_generator);
    Element facet;
    facet.type = type_facets;
    facet.element = matrix_elements.getElements(type_facets)(id);
    facet.ghost_type = gt;

    coh_model.getFacetMaterial(facet.type, facet.ghost_type)(facet.element) =
        coh_material;
    inserter.getCheckFacets(facet.type, facet.ghost_type)(facet.element) =
        false;
    insertion(facet) = true;
  }

  inserter.insertElements();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------
 */
const Array<Element>
ASRTools::insertCohElOrFacetsByCoord(const Vector<Real> & position,
                                     bool add_neighbors,
                                     bool only_double_facets) {
  AKANTU_DEBUG_IN();

  auto & mesh = model.getMesh();
  auto dim = mesh.getSpatialDimension();
  const auto gt = _not_ghost;
  const auto & pos = mesh.getNodes();
  const auto pos_it = make_view(pos, dim).begin();
  auto & node_group = mesh.createNodeGroup("doubled_nodes");

  // insert cohesive elements
  auto & coh_model = dynamic_cast<SolidMechanicsModelCohesive &>(model);
  auto & inserter = coh_model.getElementInserter();
  auto & insertion = inserter.getInsertionFacetsByElement();
  auto & mesh_facets = inserter.getMeshFacets();
  Vector<Real> bary_facet(dim);

  Real min_dist = std::numeric_limits<Real>::max();
  Element cent_facet;
  for_each_element(mesh_facets,
                   [&](auto && facet) {
                     mesh_facets.getBarycenter(facet, bary_facet);
                     auto dist = std::abs(bary_facet.distance(position));
                     if (dist < min_dist) {
                       min_dist = dist;
                       cent_facet = facet;
                     }
                   },
                   _spatial_dimension = dim - 1);

  inserter.getCheckFacets(cent_facet.type, gt)(cent_facet.element) = false;
  insertion(cent_facet) = true;
  Array<Element> facets_added;
  facets_added.push_back(cent_facet);

  if (add_neighbors) {
    auto & facet_conn = mesh_facets.getConnectivity(cent_facet.type, gt);
    CSR<Element> nodes_to_elements;
    MeshUtils::buildNode2Elements(mesh_facets, nodes_to_elements, dim - 1);

    // Synchronizer element connected to the nodes of the facet
    for (auto node : arange(2)) {
      // add corner nodes to the group
      node_group.add(facet_conn(cent_facet.element, node));

      // add corner node to the tuples
      std::get<0>(node_pairs(node)) = facet_conn(cent_facet.element, node);

      // vector of the central facet
      Vector<Real> cent_facet_dir(pos_it[facet_conn(cent_facet.element, !node)],
                                  true);
      cent_facet_dir -=
          Vector<Real>(pos_it[facet_conn(cent_facet.element, node)]);
      cent_facet_dir /= cent_facet_dir.norm();
      Real min_dot = std::numeric_limits<Real>::max();
      Element neighbor;
      Vector<Real> neighbor_facet_dir(dim);
      for (auto & elem :
           nodes_to_elements.getRow(facet_conn(cent_facet.element, node))) {
        if (elem.element == cent_facet.element)
          continue;
        if (elem.type != cent_facet.type)
          continue;
        // vector of the neighboring facet
        if (facet_conn(elem.element, 0) ==
            facet_conn(cent_facet.element, node)) {
          neighbor_facet_dir = pos_it[facet_conn(elem.element, 1)];
        } else if (facet_conn(elem.element, 1) ==
                   facet_conn(cent_facet.element, node)) {
          neighbor_facet_dir = pos_it[facet_conn(elem.element, 0)];
        } else
          AKANTU_EXCEPTION("Neighboring facet doesn't have node in common with "
                           "the central one.");
        neighbor_facet_dir -=
            Vector<Real>(pos_it[facet_conn(cent_facet.element, node)]);

        neighbor_facet_dir /= neighbor_facet_dir.norm();

        Real dot = cent_facet_dir.dot(neighbor_facet_dir);
        if (dot < min_dot) {
          min_dot = dot;
          neighbor = elem;
        }
      }

      inserter.getCheckFacets(neighbor.type, gt)(neighbor.element) = false;
      insertion(neighbor) = true;
      facets_added.push_back(neighbor);
    }
  }

  inserter.insertElements(only_double_facets);

  // flag internal facets
  auto & doubled_facets = mesh_facets.getData<bool>("doubled_facets");
  doubled_facets.initialize(mesh_facets, _spatial_dimension = dim - 1,
                            _with_nb_element = true, _default_value = false);

  for (auto & facet : facets_added) {
    Array<bool> & doubled_facets_array =
        doubled_facets(facet.type, facet.ghost_type);
    doubled_facets_array(facet.element) = true;
  }

  // create empty element group with nodes to apply Dirichlet
  model.getMesh().createElementGroupFromNodeGroup("doubled_nodes",
                                                  "doubled_nodes", dim - 1);
  AKANTU_DEBUG_OUT();
  return facets_added;
}
/* --------------------------------------------------------------------------
 */
void ASRTools::onElementsAdded(const Array<Element> & elements,
                               const NewElementsEvent &) {

  if (this->doubled_facets_ready)
    return;
  auto & mesh = model.getMesh();
  auto dim = mesh.getSpatialDimension();
  auto & mesh_facets = mesh.getMeshFacets();
  auto & doubled_facets = mesh_facets.getData<bool>("doubled_facets");
  doubled_facets.initialize(mesh_facets, _spatial_dimension = dim - 1,
                            _with_nb_element = true, _default_value = false);

  for (auto elements_range : akantu::MeshElementsByTypes(elements)) {
    auto type = elements_range.getType();
    auto ghost_type = elements_range.getGhostType();

    if (mesh.getKind(type) != _ek_regular)
      continue;

    if (ghost_type != _not_ghost)
      continue;

    auto & doubled_facets_array =
        mesh_facets.getData<bool>("doubled_facets", type, ghost_type);
    auto & element_ids = elements_range.getElements();

    for (auto && el : element_ids) {
      doubled_facets_array(el) = true;
    }
  }
  this->doubled_facets_ready = true;
}
/* --------------------------------------------------------------------------
 */
void ASRTools::onNodesAdded(const Array<UInt> & new_nodes,
                            const NewNodesEvent &) {
  AKANTU_DEBUG_IN();
  if (new_nodes.size() == 0)
    return;
  if (this->doubled_nodes_ready)
    return;
  auto & mesh = model.getMesh();
  auto & node_group = mesh.getNodeGroup("doubled_nodes");
  auto pos_it = make_view(mesh.getNodes(), mesh.getSpatialDimension()).begin();

  for (auto & node : new_nodes) {
    node_group.add(node);

    // complete pairs with lower nodes
    Vector<Real> lower_coord = pos_it[node];
    for (auto & pair : node_pairs) {
      Vector<Real> upper_coord = pos_it[std::get<0>(pair)];
      if (upper_coord == lower_coord)
        std::get<1>(pair) = node;
    }
  }

  this->doubled_nodes_ready = true;
  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
void ASRTools::applyDeltaU(Real delta_u) {
  // get element group with nodes to apply Dirichlet
  auto & crack_facets = model.getMesh().getElementGroup("doubled_nodes");

  DeltaU delta_u_bc(model, delta_u, getNodePairs());
  model.applyBC(delta_u_bc, crack_facets);
}

/* --------------------------------------------------------------------------
 */
void ASRTools::applyGelStrain(const Matrix<Real> & prestrain) {
  AKANTU_DEBUG_IN();
  auto & mesh = model.getMesh();
  auto dim = mesh.getSpatialDimension();
  AKANTU_DEBUG_ASSERT(dim == 2, "This is 2D only!");

  /// apply the new eigenstrain
  for (auto element_type :
       mesh.elementTypes(dim, _not_ghost, _ek_not_defined)) {
    Array<Real> & prestrain_vect =
        const_cast<Array<Real> &>(model.getMaterial("gel").getInternal<Real>(
            "eigen_grad_u")(element_type));
    auto prestrain_it = prestrain_vect.begin(dim, dim);
    auto prestrain_end = prestrain_vect.end(dim, dim);

    for (; prestrain_it != prestrain_end; ++prestrain_it)
      (*prestrain_it) = prestrain;
  }
  AKANTU_DEBUG_OUT();
}

} // namespace akantu
