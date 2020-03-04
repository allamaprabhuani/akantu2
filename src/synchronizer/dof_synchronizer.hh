/**
 * @file   dof_synchronizer.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 17 2011
 * @date last modification: Wed Mar 04 2020
 *
 * @brief  Synchronize Array of DOFs
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
#include "aka_array.hh"
#include "aka_common.hh"
#include "synchronizer.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
class Mesh;
class DOFManagerDefault;
} // namespace akantu

#ifndef AKANTU_DOF_SYNCHRONIZER_HH_
#define AKANTU_DOF_SYNCHRONIZER_HH_

namespace akantu {

class DOFSynchronizer : public SynchronizerImpl<Idx> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DOFSynchronizer(DOFManagerDefault & dof_manager,
                  const ID & id = "dof_synchronizer");
  ~DOFSynchronizer() override;

  virtual void registerDOFs(const ID & dof_id);
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void onNodesAdded(const Array<Idx> & nodes);

protected:
  Int getRank(const Idx & /*node*/) const final { AKANTU_TO_IMPLEMENT(); }

  /// list the entities to send to root process
  void fillEntityToSend(Array<Idx> & dofs_to_send) override;

  inline Int canScatterSize() override;
  inline Int gatheredSize() override;

  inline Idx localToGlobalEntity(const Idx & local) override;

private:
  /// information on the dofs
  DOFManagerDefault & dof_manager;
};

} // namespace akantu

#include "dof_synchronizer_inline_impl.hh"

#endif /* AKANTU_DOF_SYNCHRONIZER_HH_ */
