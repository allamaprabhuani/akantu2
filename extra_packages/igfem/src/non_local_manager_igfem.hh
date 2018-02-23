/**
 * @file   non_local_manager_igfem.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Sep 21 14:21:33 2015
 *
 * @brief Class that manages all the non-local neighborhoods for IGFEM
 * simulations
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
#ifdef AKANTU_DAMAGE_NON_LOCAL
#ifndef __AKANTU_NON_LOCAL_MANAGER_IGFEM_HH__
#define __AKANTU_NON_LOCAL_MANAGER_IGFEM_HH__
/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model_igfem.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class NonLocalManagerIGFEM : public NonLocalManager {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLocalManagerIGFEM(SolidMechanicsModelIGFEM & model,
                       const ID & id = "non_local_manager_igfem",
                       const MemoryID & memory_id = 0);
  virtual ~NonLocalManagerIGFEM();

  /* --------------------------------------------------------------------------
   */
  /* Methods */
  /* --------------------------------------------------------------------------
   */
public:
  /// initialize the non-local manager: compute pair lists and weights for all
  /// neighborhoods
  virtual void init();

  /// average the internals and compute the non-local stresses
  virtual void computeAllNonLocalStresses();

  /* --------------------------------------------------------------------------
   */
  /* MeshEventHandler inherited members */
  /* --------------------------------------------------------------------------
   */

  virtual void
  onElementsRemoved(const Array<Element> & element_list,
                    const ElementTypeMapArray<UInt> & new_numbering,
                    const RemovedElementsEvent & event);

  virtual void onElementsAdded(const Array<Element> & element_list,
                               const NewElementsEvent & event);

private:
  /// cleanup unneccessary ghosts
  virtual void
  cleanupExtraGhostElements(ElementTypeMap<UInt> & nb_ghost_protected);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
};

__END_AKANTU__

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#endif /* __AKANTU_NON_LOCAL_MANAGER_IGFEM_HH__ */
#endif /* AKANTU_DAMAGE_NON_LOCAL */
