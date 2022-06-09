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
  void findContactNodes();

  // TODO: move to a new class
  Matrix<Real> constructInterpolationMatrix(NodeGroup & ref_node_group,
                                            NodeGroup & eval_node_group,
                                            Array<Real> eval_radiuses);
  /// construct Phi matrix used for interpolation
  Matrix<Real> constructPhiMatrix(NodeGroup & ref_node_group,
                                  NodeGroup & eval_node_group,
                                  Array<Real> & eval_radiuses);

private:
  /// reads the input file to get contact detection options
  void parseSection(const ParserSection & section) override;

  /// compute radius to detect contact nodes
  std::pair<Array<Real>, Array<UInt>>
  computeRadiuses(NodeGroup & ref_node_group, NodeGroup & eval_node_group);

  /// radial basis function
  Real computeRadialBasisInterpolation(Real distance, Real radius);

  /// distances between a reference node and a other nodes
  Array<Real> computeDistancesToRefNode(UInt & ref_node,
                                        NodeGroup & eval_node_group);

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
