/**
 * @file   contact_detector_internodes.hh
 *
 * @author Moritz Waldleben <moritz.waldleben@epfl.ch>
 *
 * @date creation: Thu Jul 09 2022
 * @date last modification: Wed Dec 14 2022
 *
 * @brief Algorithm to detetect contact nodes
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "parsable.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_DETECTOR_INTERNODES_HH__
#define __AKANTU_CONTACT_DETECTOR_INTERNODES_HH__

namespace akantu {
class Mesh;
class NodeGroup;
} // namespace akantu

namespace akantu {

/* -------------------------------------------------------------------------- */
class ContactDetectorInternodes : public Parsable {

  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructors                                                  */
  /* ------------------------------------------------------------------------ */
public:
  ContactDetectorInternodes(Mesh & mesh, const ID & id = "contact_detector");

  ~ContactDetectorInternodes() override = default;

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// find contact nodes for iteration
  void findContactNodes(Real C = 0.95);

  /// construct interpolation matrices
  Matrix<Real> constructInterpolationMatrix(const NodeGroup & ref_node_group,
                                            const NodeGroup & eval_node_group,
                                            Array<Real> eval_radiuses);

  /// construct Phi matrix used for interpolation
  Matrix<Real> constructPhiMatrix(const NodeGroup & ref_node_group,
                                  const NodeGroup & eval_node_group,
                                  Array<Real> & eval_radiuses);

private:
  /// compute radius to detect contact nodes
  void computeRadiuses(Array<Real> & radius_parameters,
                       const NodeGroup & node_group, Real c = 0.5,
                       Real C = 0.95);

  /// reads the input file to get contact detection options
  void parseSection(const ParserSection & section) override;

  /// radial basis function
  Real evaluateRadialBasisFunction(Real distance, Real radius);

  /// distances between a reference node and other nodes
  /// the out_array must already be allocated with sufficient size
  void computeDistancesToRefNode(UInt & ref_node,
                                 const NodeGroup & eval_node_group,
                                 Array<Real> & out_array);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the mesh
  AKANTU_GET_MACRO(Mesh, mesh, Mesh &)

  /// get radiuses of attack of master nodes
  AKANTU_GET_MACRO_NOT_CONST(MasterRadiuses, master_radiuses, Array<Real> &)

  /// get radiuses of attack of master nodes
  AKANTU_GET_MACRO_NOT_CONST(SlaveRadiuses, slave_radiuses, Array<Real> &)

  /// get initial master node group
  NodeGroup & getInitialMasterNodeGroup();

  /// get initial slave node group
  NodeGroup & getInitialSlaveNodeGroup();

  /// get master node group
  NodeGroup & getMasterNodeGroup();

  /// get slave node group
  NodeGroup & getSlaveNodeGroup();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// mesh
  Mesh & mesh;

  /// dimension of the model
  UInt spatial_dimension{0};

  /// id of master node group, i.e nodes on interface
  ID id_master_nodes{};

  /// id of slave nodes group, i.e nodes on interface
  ID id_slave_nodes{};

  /// attack radiuses for master nodes
  Array<Real> master_radiuses{0};

  /// attack radiuses for master nodes
  Array<Real> slave_radiuses{0};

  /// blocked boundary dofs array
  std::unique_ptr<Array<Real>> blocked_dofs;
};

} // namespace akantu

#endif /* __AKANTU_CONTACT_DETECTOR_INTERNODES_HH__ */
