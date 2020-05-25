/**
 * @file   solid_mechanics_model_RVE.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Jan 13 15:32:35 2016
 *
 * @brief  Implementation of SolidMechanicsModelRVE
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
#include "solid_mechanics_model_RVE.hh"
#include "aka_random_generator.hh"
#include "element_group.hh"
#include "material_damage_iterative.hh"
#include "node_group.hh"
#include "non_linear_solver.hh"
#include "non_local_manager.hh"
#include "parser.hh"
#include "sparse_matrix.hh"
#include <string>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SolidMechanicsModelRVE::SolidMechanicsModelRVE(Mesh & mesh,
                                               bool use_RVE_mat_selector,
                                               UInt nb_gel_pockets, UInt dim,
                                               const ID & id,
                                               const MemoryID & memory_id)
    : SolidMechanicsModel(mesh, dim, id, memory_id),
      ASRTools(dynamic_cast<SolidMechanicsModel &>(*this)),
      use_RVE_mat_selector(use_RVE_mat_selector),
      nb_gel_pockets(nb_gel_pockets), stiffness_changed(true) {
  AKANTU_DEBUG_IN();
  RandomGenerator<UInt>::seed(1);

  /// remove the corner nodes from the surface node groups:
  /// This most be done because corner nodes a not periodic
  mesh.getElementGroup("top").removeNode(corner_nodes(2));
  mesh.getElementGroup("top").removeNode(corner_nodes(3));
  mesh.getElementGroup("left").removeNode(corner_nodes(3));
  mesh.getElementGroup("left").removeNode(corner_nodes(0));
  mesh.getElementGroup("bottom").removeNode(corner_nodes(1));
  mesh.getElementGroup("bottom").removeNode(corner_nodes(0));
  mesh.getElementGroup("right").removeNode(corner_nodes(2));
  mesh.getElementGroup("right").removeNode(corner_nodes(1));

  const auto & bottom = mesh.getElementGroup("bottom").getNodeGroup();
  bottom_nodes.insert(bottom.begin(), bottom.end());

  const auto & left = mesh.getElementGroup("left").getNodeGroup();
  left_nodes.insert(left.begin(), left.end());

  mesh.makePeriodic(_x, "left", "right");
  mesh.makePeriodic(_y, "bottom", "top");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolidMechanicsModelRVE::~SolidMechanicsModelRVE() = default;

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::initFullImpl(const ModelOptions & options) {
  AKANTU_DEBUG_IN();

  auto options_cp(options);
  options_cp.analysis_method = AnalysisMethod::_static;

  SolidMechanicsModel::initFullImpl(options_cp);

  auto & solver = this->getNonLinearSolver();
  solver.set("max_iterations", 50);
  solver.set("threshold", 1e-5);
  solver.set("convergence_type", SolveConvergenceCriteria::_solution);
  // this->initMaterials();

  /// variables for parallel execution
  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();

  /// compute the volume of the RVE
  this->computeModelVolume();
  if (prank == 0)
    std::cout << "The volume of the RVE is " << this->volume << std::endl;

  /// dumping
  std::stringstream base_name;
  base_name << this->id; // << this->memory_id - 1;
  this->setBaseName(base_name.str());
  this->addDumpFieldVector("displacement");
  this->addDumpField("stress");
  this->addDumpField("grad_u");
  this->addDumpField("eigen_grad_u");
  this->addDumpField("blocked_dofs");
  this->addDumpField("material_index");
  this->addDumpField("damage");
  this->addDumpField("Sc");
  this->addDumpField("external_force");
  this->addDumpField("equivalent_stress");
  this->addDumpField("internal_force");
  this->addDumpField("delta_T");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::assembleInternalForces() {

  /// displacement correction
  auto disp_begin = make_view(*this->displacement, spatial_dimension).begin();

  auto u_top_corr = Vector<Real>(disp_begin[this->corner_nodes(2)]) -
                    Vector<Real>(disp_begin[this->corner_nodes(1)]);
  auto u_right_corr = Vector<Real>(disp_begin[this->corner_nodes(2)]) -
                      Vector<Real>(disp_begin[this->corner_nodes(3)]);

  auto correct_disp = [&](auto && node, auto && correction) {
    const auto & pair = mesh.getPeriodicMasterSlaves().equal_range(node);
    for (auto && data : range(pair.first, pair.second)) {
      auto slave = data.second;

      auto slave_disp = Vector<Real>(disp_begin[slave]);
      slave_disp = Vector<Real>(disp_begin[node]);
      slave_disp += correction;
    }
  };

  for (auto node : bottom_nodes) {
    correct_disp(node, u_top_corr);
  }

  for (auto node : left_nodes) {
    correct_disp(node, u_right_corr);
  }

  SolidMechanicsModel::assembleInternalForces();
}

/* --------------------------------------------------------------------------
 */
void SolidMechanicsModelRVE::advanceASR(const Matrix<Real> & prestrain) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(spatial_dimension == 2, "This is 2D only!");

  /// apply the new eigenstrain
  for (auto element_type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_not_defined)) {
    Array<Real> & prestrain_vect =
        const_cast<Array<Real> &>(this->getMaterial("gel").getInternal<Real>(
            "eigen_grad_u")(element_type));
    auto prestrain_it =
        prestrain_vect.begin(spatial_dimension, spatial_dimension);
    auto prestrain_end =
        prestrain_vect.end(spatial_dimension, spatial_dimension);

    for (; prestrain_it != prestrain_end; ++prestrain_it)
      (*prestrain_it) = prestrain;
  }

  /// advance the damage
  MaterialDamageIterativeInterface & mat_paste =
      dynamic_cast<MaterialDamageIterativeInterface &>(
          this->getMaterial("paste"));
  MaterialDamageIterativeInterface & mat_aggregate =
      dynamic_cast<MaterialDamageIterativeInterface &>(
          this->getMaterial("aggregate"));
  UInt nb_damaged_elements = 0;
  Real max_eq_stress_aggregate = 0;
  Real max_eq_stress_paste = 0;
  this->stiffness_changed = false;

  /// variables for parallel execution
  auto && comm = akantu::Communicator::getWorldCommunicator();
  auto prank = comm.whoAmI();

  /// save nodal fields (disp, boun & ext_force)
  this->storeNodalFields();

  do {
    /// restore nodals and update grad_u accordingly
    this->restoreNodalFields();

    /// restore historical internal fields (sigma_v for visc)
    this->restoreInternalFields();

    this->solveStep();

    /// compute damage
    max_eq_stress_aggregate = mat_aggregate.getNormMaxEquivalentStress();
    max_eq_stress_paste = mat_paste.getNormMaxEquivalentStress();

    nb_damaged_elements = 0;
    if (max_eq_stress_aggregate > max_eq_stress_paste)
      nb_damaged_elements = mat_aggregate.updateDamage();
    else if (max_eq_stress_aggregate < max_eq_stress_paste)
      nb_damaged_elements = mat_paste.updateDamage();
    else
      nb_damaged_elements =
          (mat_paste.updateDamage() + mat_aggregate.updateDamage());

    /// mark the flag to update stiffness if elements were damaged
    if (nb_damaged_elements)
      this->stiffness_changed = true;

    std::cout << "Proc " << prank << " the number of damaged elements is "
              << nb_damaged_elements << std::endl;

  } while (nb_damaged_elements);

  //  if (this->nb_dumps % 10 == 0) {
  // this->dump();
  //  }
  // this->nb_dumps += 1;

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
 */
void SolidMechanicsModelRVE::initMaterials() {
  AKANTU_DEBUG_IN();

  // make sure the material are instantiated
  if (!are_materials_instantiated)
    instantiateMaterials();

  if (use_RVE_mat_selector) {
    const Vector<Real> & lowerBounds = mesh.getLowerBounds();
    const Vector<Real> & upperBounds = mesh.getUpperBounds();
    Real bottom = lowerBounds(1);
    Real top = upperBounds(1);
    Real box_size = std::abs(top - bottom);
    Real eps = box_size * 1e-6;

    auto tmp = std::make_shared<GelMaterialSelector>(
        *this, box_size, "gel", this->nb_gel_pockets, "aggregate", eps);
    tmp->setFallback(material_selector);
    material_selector = tmp;
  }

  this->assignMaterialToElements();
  // synchronize the element material arrays
  this->synchronize(SynchronizationTag::_material_id);

  for (auto & material : materials) {
    /// init internals properties
    const auto type = material->getID();
    if (type.find("material_FE2") != std::string::npos)
      continue;
    material->initMaterial();
  }

  this->synchronize(SynchronizationTag::_smm_init_mat);

  if (this->non_local_manager) {
    this->non_local_manager->initialize();
  }
  // SolidMechanicsModel::initMaterials();

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
