/**
 * @file test_embedded_interface_model_prestress.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: mar. avril 28 2015
 * @date last modification: mar. avril 28 2015
 *
 * @brief Embedded model test for prestressing (bases on stress norm)
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "aka_common.hh"
#include "embedded_interface_model.hh"

/* -------------------------------------------------------------------------- */

using namespace akantu;

#define YG 0.483644859
#define I_eq 0.012488874
#define A_eq (1e-2 + 1. / 7. * 1.)

/* -------------------------------------------------------------------------- */

struct StressSolution : public BC::Neumann::FromHigherDim {
  Real M;
  Real I;
  Real yg;
  Real pre_stress;

  StressSolution(UInt dim, Real M, Real I, Real yg = 0, Real pre_stress = 0) :
    BC::Neumann::FromHigherDim(Matrix<Real>(dim, dim)),
    M(M), I(I), yg(yg), pre_stress(pre_stress)
  {}

  void operator()(const QuadraturePoint & quad_point,
                          Vector<Real> & dual,
                          const Vector<Real> & coord,
                          const Vector<Real> & normals) const {
    UInt dim = coord.size();

    if (dim < 2) AKANTU_DEBUG_ERROR("Solution not valid for 1D");

    Matrix<Real> stress(dim, dim); stress.clear();
    stress(0, 0) = this->stress(coord(1));
    dual.mul<false>(stress, normals);
  }

  inline Real stress(Real height) const {
    return -M / I * (height - yg) + pre_stress;
  }

  inline Real neutral_axis() const {
    return -I * pre_stress / M + yg;
  }
};

/* -------------------------------------------------------------------------- */

int main (int argc, char * argv[]) {
  debug::setDebugLevel(dblWarning);
  initialize("prestress.dat", argc, argv);

  Math::setTolerance(1e-6);

  const UInt dim = 2;

/* -------------------------------------------------------------------------- */

  Mesh mesh(dim);
  mesh.read("embedded_mesh_prestress.msh");
  mesh.createGroupsFromMeshData<std::string>("physical_names");

  EmbeddedInterfaceModel model(mesh, dim);
  model.initFull(SolidMechanicsModelOptions(_static));

/* -------------------------------------------------------------------------- */
/* Computation of analytical residual                                         */
/* -------------------------------------------------------------------------- */

  /*
   * q = 1000 N/m
   * L = 20 m
   * a = 1 m
   */

  Real steel_area = model.getMaterial("elastic_r").getParam<Real>("area");
  Real pre_stress = 1e6;
  Real stress_norm = 0.;

  StressSolution * concrete_stress = NULL, * steel_stress = NULL;

  Real pre_force = pre_stress * steel_area;
  Real pre_moment = -pre_force * (YG - 0.25);
  Real neutral_axis = YG - I_eq / A_eq * pre_force / pre_moment;

  concrete_stress = new StressSolution(dim, pre_moment, 7. * I_eq, YG, -pre_force / (7. * A_eq));
  steel_stress = new StressSolution(dim, pre_moment, I_eq, YG, pre_stress - pre_force / A_eq);

  stress_norm = std::abs(concrete_stress->stress(1)) * (1 - neutral_axis) * 0.5
    + std::abs(concrete_stress->stress(0)) * neutral_axis * 0.5
    + std::abs(steel_stress->stress(0.25)) * steel_area;

  model.applyBC(*concrete_stress, "XBlocked");
  ElementGroup & end_node = mesh.getElementGroup("EndNode");
  NodeGroup & end_node_group = end_node.getNodeGroup();
  NodeGroup::const_node_iterator end_node_it = end_node_group.begin();

  Vector<Real> end_node_force = model.getForce().begin(dim)[*end_node_it];
  end_node_force(0) += steel_stress->stress(0.25) * steel_area; 

  Array<Real> analytical_residual(mesh.getNbNodes(), dim, "analytical_residual");
  analytical_residual.copy(model.getForce());
  model.getForce().clear();

  delete concrete_stress;
  delete steel_stress;

/* -------------------------------------------------------------------------- */

  model.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "XBlocked");
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "YBlocked");

  // Assemble the global stiffness matrix
  model.assembleStiffnessMatrix();
  
  model.updateResidual();

  if (!model.solveStatic<_scm_newton_raphson_tangent_not_computed, _scc_residual>(1e-6, 1))
    return EXIT_FAILURE;

  model.updateResidual();

/* -------------------------------------------------------------------------- */
/* Computation of FEM residual norm                                           */
/* -------------------------------------------------------------------------- */

  ElementGroup & xblocked = mesh.getElementGroup("XBlocked");
  NodeGroup & boundary_nodes = xblocked.getNodeGroup();
  NodeGroup::const_node_iterator
    nodes_it = boundary_nodes.begin(),
    nodes_end = boundary_nodes.end();
  Array<Real>::vector_iterator com_res = model.getResidual().begin(dim);
  Array<Real>::vector_iterator ana_res = analytical_residual.begin(dim);
  Array<Real>::vector_iterator position = mesh.getNodes().begin(dim);

  Real res_sum = 0.;
  UInt lower_node = -1;
  UInt upper_node = -1;
  Real lower_dist = 1;
  Real upper_dist = 1;

  for (; nodes_it != nodes_end ; ++nodes_it) {
    UInt node_number = *nodes_it;
    const Vector<Real> res = com_res[node_number];
    const Vector<Real> ana = ana_res[node_number];
    const Vector<Real> pos = position[node_number];

    if (!Math::are_float_equal(pos(1), 0.25)) {
      if ((std::abs(pos(1) - 0.25) < lower_dist) && (pos(1) < 0.25)) {
        lower_dist = std::abs(pos(1) - 0.25);
        lower_node = node_number;
      }

      if ((std::abs(pos(1) - 0.25) < upper_dist) && (pos(1) > 0.25)) {
        upper_dist = std::abs(pos(1) - 0.25);
        upper_node = node_number;
      }
    }

    for (UInt i = 0 ; i < dim ; i++) {
      if (!Math::are_float_equal(pos(1), 0.25)) {
        res_sum += std::abs(res(i));
      }
    }
  }

  const Vector<Real> upper_res = com_res[upper_node], lower_res = com_res[lower_node];
  const Vector<Real> end_node_res = com_res[*end_node_it];
  Vector<Real> delta = upper_res - lower_res;
  delta *= lower_dist / (upper_dist + lower_dist);
  Vector<Real> concrete_residual = lower_res + delta;
  Vector<Real> steel_residual = end_node_res - concrete_residual;

  for (UInt i = 0 ; i < dim ; i++) {
    res_sum += std::abs(concrete_residual(i));
    res_sum += std::abs(steel_residual(i));
  }

  if (std::abs(res_sum - stress_norm) / stress_norm > 1e-3)
    return EXIT_FAILURE;

  finalize();
  return 0;
}
