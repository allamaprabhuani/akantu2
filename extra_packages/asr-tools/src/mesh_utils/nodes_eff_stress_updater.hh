/**
 * @file   nodes_eff_stress_updater.hh
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @date Tue Feb 8  2022
 * @brief  synchronizer for the effective stresses at the shared nodes
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
#ifndef __AKANTU_NODES_EFF_STRESS_UPDATER_HH__
#define __AKANTU_NODES_EFF_STRESS_UPDATER_HH__

/* -------------------------------------------------------------------------- */
#include "data_accessor.hh"
#include "mesh.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {
class ElementSynchronizer;
}

namespace akantu {

class NodesEffStressUpdater : public DataAccessor<Element> {
public:
  NodesEffStressUpdater(Mesh & mesh, ElementSynchronizer & synchronizer,
                        Array<Real> & nodes_eff_stress)
      : mesh(mesh), synchronizer(synchronizer),
        nodes_eff_stress(nodes_eff_stress) {
    AKANTU_DEBUG_ASSERT(mesh.getNbNodes() == nodes_eff_stress.size(),
                        "Array nodes_eff_stress does not have the same size as "
                        "the number of nodes in the mesh");
  }

  void updateMaxEffStressAtNodes();

  /* ----------------------------------------------------------------- */
  /* Data Accessor inherited members                                   */
  /* ----------------------------------------------------------------- */
public:
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  /* ----------------------------------------------------------------- */
  /* Members                                                           */
  /* ----------------------------------------------------------------- */
private:
  /// Reference to the mesh to update
  Mesh & mesh;

  /// distributed synchronizer to communicate nodes positions
  ElementSynchronizer & synchronizer;

  /// average effective stress of facets loop around a node
  Array<Real> & nodes_eff_stress;
};

} // namespace akantu

#include "nodes_eff_stress_updater_inline_impl.hh"

#endif /* __AKANTU_NODES_EFF_STRESS_UPDATER_HH__ */
