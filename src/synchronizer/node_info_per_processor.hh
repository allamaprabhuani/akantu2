/**
 * Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "communication_buffer.hh"
#include "mesh_accessor.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NODE_INFO_PER_PROCESSOR_HH_
#define AKANTU_NODE_INFO_PER_PROCESSOR_HH_

namespace akantu {
class NodeSynchronizer;
class Communicator;
} // namespace akantu

/* -------------------------------------------------------------------------- */

namespace akantu {

class NodeInfoPerProc : protected MeshAccessor {
public:
  NodeInfoPerProc(NodeSynchronizer & synchronizer, Int message_cnt, Int root);

  void synchronize();

protected:
  virtual void synchronizeNodes() = 0;
  virtual void synchronizeTypes() = 0;
  virtual void synchronizeGroups() = 0;
  virtual void synchronizePeriodicity() = 0;
  virtual void synchronizeTags() = 0;

protected:
  template <class CommunicationBuffer>
  void fillNodeGroupsFromBuffer(CommunicationBuffer & buffer);
  void fillNodesType();

  void fillCommunicationScheme(const Array<Idx> &);
  void fillNodalData(DynamicCommunicationBuffer & buffer,
                     const std::string & tag_name);

  void fillPeriodicPairs(const Array<Idx> &, std::vector<Idx> &);
  void receiveMissingPeriodic(DynamicCommunicationBuffer &);

protected:
  NodeSynchronizer & synchronizer;
  const Communicator & comm;
  Int rank;
  Int nb_proc;
  Int root;

  Mesh & mesh;

  Int spatial_dimension;
  Int message_count;
};

/* -------------------------------------------------------------------------- */
class MasterNodeInfoPerProc : public NodeInfoPerProc {
public:
  MasterNodeInfoPerProc(NodeSynchronizer & synchronizer, Int message_cnt,
                        Int root);

  void synchronizeNodes() override;
  void synchronizeTypes() override;
  void synchronizeGroups() override;
  void synchronizePeriodicity() override;
  void synchronizeTags() override;

private:
  void fillTagBuffers(std::vector<DynamicCommunicationBuffer> & buffers,
                      const std::string & tag_name);

  /// get the list of nodes to send and send them
  std::vector<Array<Idx>> nodes_per_proc;
  Array<Int> nb_nodes_per_proc;
  Array<Real> all_nodes;
  Array<NodeFlag> all_periodic_flags;
  Array<Int> nodes_pranks;
};

/* -------------------------------------------------------------------------- */
class SlaveNodeInfoPerProc : public NodeInfoPerProc {
public:
  SlaveNodeInfoPerProc(NodeSynchronizer & synchronizer, Int message_cnt,
                       Int root);

  void synchronizeNodes() override;
  void synchronizeTypes() override;
  void synchronizeGroups() override;
  void synchronizePeriodicity() override;
  void synchronizeTags() override;

private:
};

} // namespace akantu

#endif /* AKANTU_NODE_INFO_PER_PROCESSOR_HH_ */
