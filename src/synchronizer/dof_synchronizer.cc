/**
 * @file   dof_synchronizer.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 17 2011
 * @date last modification: Wed Oct 21 2015
 *
 * @brief  DOF synchronizing object implementation
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "dof_synchronizer.hh"
#include "aka_zip.hh"
#include "dof_manager_default.hh"
#include "mesh.hh"
#include "node_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
/**
 * A DOFSynchronizer needs a mesh and the number of degrees of freedom
 * per node to be created. In the constructor computes the local and global dof
 * number for each dof. The member
 * proc_informations (std vector) is resized with the number of mpi
 * processes. Each entry in the vector is a PerProcInformations object
 * that contains the interactions of the current mpi process (prank) with the
 * mpi process corresponding to the position of that entry. Every
 * ProcInformations object contains one array with the dofs that have
 * to be sent to prank and a second one with dofs that willl be received form
 * prank.
 * This information is needed for the asychronous communications. The
 * constructor sets up this information.
 */
DOFSynchronizer::DOFSynchronizer(DOFManagerDefault & dof_manager, const ID & id,
                                 MemoryID memory_id,
                                 const StaticCommunicator & comm)
    : SynchronizerImpl<UInt>(id, memory_id, comm), root(0),
      dof_manager(dof_manager), root_dofs(0, 1, "dofs-to-receive-from-master"),
      dof_changed(true) {
  std::vector<ID> dof_ids = dof_manager.getDOFIDs();

  // Transfers nodes to global equation numbers in new schemes
  for (ID dof_id : dof_ids) {
    registerDOFs(dof_id);
  }

  this->initScatterGatherCommunicationScheme();
}

/* -------------------------------------------------------------------------- */
DOFSynchronizer::~DOFSynchronizer() {}

/* -------------------------------------------------------------------------- */
void DOFSynchronizer::registerDOFs(const ID & dof_id) {
  if (this->nb_proc == 1)
    return;

  if (dof_manager.getSupportType(dof_id) != _dst_nodal)
    return;

  using const_scheme_iterator = Communications<UInt>::const_scheme_iterator;

  const auto equation_numbers = dof_manager.getLocalEquationNumbers(dof_id);

  const auto & associated_nodes = dof_manager.getDOFsAssociatedNodes(dof_id);
  const auto & node_synchronizer = dof_manager.getMesh().getNodeSynchronizer();
  const auto & node_communications = node_synchronizer.getCommunications();

  auto transcode_node_to_global_dof_scheme =
      [this, &associated_nodes,
       &equation_numbers](const_scheme_iterator it, const_scheme_iterator end,
                          const CommunicationSendRecv & sr) -> void {
    for (; it != end; ++it) {
      auto & scheme = communications.createScheme(it->first, sr);

      const auto & node_scheme = it->second;
      for (auto & node : node_scheme) {
        auto an_begin = associated_nodes.begin();
        auto an_it = an_begin;
        auto an_end = associated_nodes.end();

        std::vector<UInt> global_dofs_per_node;
        while ((an_it = std::find(an_it, an_end, node)) != an_end) {
          UInt pos = an_it - an_begin;
          UInt local_eq_num = equation_numbers(pos);
          UInt global_eq_num =
              dof_manager.localToGlobalEquationNumber(local_eq_num);
          global_dofs_per_node.push_back(global_eq_num);
          ++an_it;
        }

        std::sort(global_dofs_per_node.begin(), global_dofs_per_node.end());
        std::transform(global_dofs_per_node.begin(), global_dofs_per_node.end(),
                       global_dofs_per_node.begin(), [this](UInt g) -> UInt {
                         UInt l = dof_manager.globalToLocalEquationNumber(g);
                         return l;
                       });
        for (auto & leqnum : global_dofs_per_node) {
          scheme.push_back(leqnum);
        }
      }
    }
  };

  for (auto sr_it = send_recv_t::begin(); sr_it != send_recv_t::end();
       ++sr_it) {
    auto ncs_it = node_communications.begin_scheme(*sr_it);
    auto ncs_end = node_communications.end_scheme(*sr_it);

    transcode_node_to_global_dof_scheme(ncs_it, ncs_end, *sr_it);
  }

  dof_changed = true;
}

/* -------------------------------------------------------------------------- */
void DOFSynchronizer::initScatterGatherCommunicationScheme() {
  AKANTU_DEBUG_IN();

  if (this->nb_proc == 1) {
    dof_changed = false;
    AKANTU_DEBUG_OUT();
    return;
  }

  UInt nb_dofs = dof_manager.getLocalSystemSize();

  this->root_dofs.clear();
  this->master_receive_dofs.clear();

  Array<UInt> dofs_to_send;
  for (UInt n = 0; n < nb_dofs; ++n) {
    if (dof_manager.isLocalOrMasterDOF(n)) {
      auto & receive_per_proc = master_receive_dofs[this->root];
      UInt global_dof = dof_manager.localToGlobalEquationNumber(n);

      root_dofs.push_back(n);
      receive_per_proc.push_back(global_dof);
      dofs_to_send.push_back(global_dof);
    }
  }

  if (this->rank == UInt(this->root)) {
    Array<UInt> nb_dof_per_proc(this->nb_proc);
    communicator.gather(dofs_to_send.getSize(), nb_dof_per_proc);

    std::vector<CommunicationRequest> requests;
    for (UInt p = 0; p < nb_proc; ++p) {
      if (p == UInt(this->root))
        continue;

      auto & receive_per_proc = master_receive_dofs[p];
      receive_per_proc.resize(nb_dof_per_proc(p));
      if (nb_dof_per_proc(p) != 0)
        requests.push_back(communicator.asyncReceive(
            receive_per_proc, p,
            Tag::genTag(p, 0, Tag::_GATHER_INITIALIZATION, this->hash_id)));
    }

    communicator.waitAll(requests);
    communicator.freeCommunicationRequest(requests);
  } else {
    communicator.gather(dofs_to_send.getSize(), this->root);
    AKANTU_DEBUG(dblDebug, "I have " << nb_dofs << " dofs ("
                                     << dofs_to_send.getSize()
                                     << " to send to master proc");

    if (dofs_to_send.getSize() != 0)
      communicator.send(dofs_to_send, this->root,
                        Tag::genTag(this->rank, 0, Tag::_GATHER_INITIALIZATION,
                                    this->hash_id));
  }

  dof_changed = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool DOFSynchronizer::hasChanged() {
  communicator.allReduce(dof_changed, _so_lor);
  return dof_changed;
}

/* -------------------------------------------------------------------------- */
void DOFSynchronizer::onNodesAdded(const Array<UInt> & nodes_list) {
  auto dof_ids = dof_manager.getDOFIDs();

  const auto & node_synchronizer = dof_manager.getMesh().getNodeSynchronizer();
  const auto & node_communications = node_synchronizer.getCommunications();

  std::set<UInt> relevant_nodes;
  std::map<UInt, std::vector<UInt>> nodes_per_proc[2];
  for (auto sr_it = send_recv_t::begin(); sr_it != send_recv_t::end();
       ++sr_it) {
    auto sit = node_communications.begin_scheme(*sr_it);
    auto send = node_communications.end_scheme(*sr_it);

    for (; sit != send; ++sit) {
      auto proc = sit->first;
      const auto & scheme = sit->second;
      for (auto node : nodes_list) {
        if (scheme.find(node) == -1)
          continue;
        relevant_nodes.insert(node);
        nodes_per_proc[*sr_it][proc].push_back(node);
      }
    }
  }

  std::map<UInt, std::vector<UInt>> dofs_per_proc[2];
  for (auto & dof_id : dof_ids) {
    const auto & associated_nodes = dof_manager.getDOFsAssociatedNodes(dof_id);
    const auto & local_equation_numbers =
        dof_manager.getEquationsNumbers(dof_id);

    for (auto tuple : zip(associated_nodes, local_equation_numbers)) {
      UInt assoc_node;
      UInt local_eq_num;
      std::tie(assoc_node, local_eq_num) = tuple;

      for (auto sr_it = send_recv_t::begin(); sr_it != send_recv_t::end();
           ++sr_it) {
        for (auto & pair : nodes_per_proc[*sr_it]) {
          if (std::find(pair.second.end(), pair.second.end(), assoc_node) !=
              pair.second.end()) {
            dofs_per_proc[*sr_it][pair.first].push_back(local_eq_num);
          }
        }
      }
    }
  }

  for (auto sr_it = send_recv_t::begin(); sr_it != send_recv_t::end();
           ++sr_it) {
    for (auto & pair : dofs_per_proc[*sr_it]) {
      std::sort(pair.second.begin(), pair.second.end(),
                [this] (UInt la, UInt lb) -> bool {
                  UInt ga = dof_manager.localToGlobalEquationNumber(la);
                  UInt gb = dof_manager.localToGlobalEquationNumber(lb);
                  return ga < gb;
                });

      auto & scheme = communications.getScheme(pair.first, *sr_it);
      for(auto leq : pair.second) {
        scheme.push_back(leq);
      }
    }
  }
  dof_changed = true;
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
