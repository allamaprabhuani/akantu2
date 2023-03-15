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
#include "aka_common.hh"
#include "communication_buffer.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ELEMENT_INFO_PER_PROCESSOR_HH_
#define AKANTU_ELEMENT_INFO_PER_PROCESSOR_HH_

namespace akantu {
class ElementSynchronizer;
class Communicator;
class MeshPartition;
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {

class ElementInfoPerProc : protected MeshAccessor {
public:
  ElementInfoPerProc(ElementSynchronizer & synchronizer, Int message_cnt,
                     Int root, ElementType type);
  bool synchronize();

protected:
  virtual void synchronizeConnectivities() = 0;
  virtual void synchronizePartitions() = 0;
  virtual void synchronizeTags() = 0;
  virtual void synchronizeGroups() = 0;
  virtual bool needSynchronize() = 0;

protected:
  void fillCommunicationScheme(const Array<Int> & partition);

  template <class CommunicationBuffer>
  void fillElementGroupsFromBuffer(CommunicationBuffer & buffer);

  template <typename T, typename BufferType>
  void fillMeshDataTemplated(BufferType & buffer, const std::string & tag_name,
                             Int nb_component);

  template <typename BufferType>
  void fillMeshData(BufferType & buffer, const std::string & tag_name,
                    const MeshDataTypeCode & type_code, Int nb_component);

protected:
  ElementSynchronizer & synchronizer;

  Int rank{0};
  Int nb_proc{1};

  Int root{0};

  ElementType type{_not_defined};

  Int nb_tags{0};
  Int nb_nodes_per_element{0};
  Int nb_element{0};

  Int nb_local_element{0};
  Int nb_ghost_element{0};

  Int message_count{0};
  Mesh & mesh;
  const Communicator & comm;
};

/* -------------------------------------------------------------------------- */
class MasterElementInfoPerProc : public ElementInfoPerProc {
public:
  MasterElementInfoPerProc(ElementSynchronizer & synchronizer, Int message_cnt,
                           Int root, ElementType type,
                           const MeshPartition & partition);

protected:
  void synchronizeConnectivities() override;
  void synchronizePartitions() override;
  void synchronizeTags() override;
  void synchronizeGroups() override;
  bool needSynchronize() override { return type != _not_defined; }

protected:
  template <typename T>
  void fillTagBufferTemplated(std::vector<DynamicCommunicationBuffer> & buffers,
                              const std::string & tag_name);
  void fillTagBuffer(std::vector<DynamicCommunicationBuffer> & buffers,
                     const std::string & tag_name);

private:
  const MeshPartition & partition;

  Vector<Int> all_nb_local_element;
  Vector<Int> all_nb_ghost_element;
  Vector<Int> all_nb_element_to_send;
};

/* -------------------------------------------------------------------------- */
class SlaveElementInfoPerProc : public ElementInfoPerProc {
public:
  SlaveElementInfoPerProc(ElementSynchronizer & synchronizer, Int message_cnt,
                          Int root);

protected:
  void synchronizeConnectivities() override;
  void synchronizePartitions() override;
  void synchronizeTags() override;
  void synchronizeGroups() override;

  bool needSynchronize() override;

private:
  Int nb_element_to_receive{0};
};

} // namespace akantu

#include "element_info_per_processor_tmpl.hh"

#endif /* AKANTU_ELEMENT_INFO_PER_PROCESSOR_HH_ */
