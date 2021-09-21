/**
 * @file   facet_synchronizer.hh
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Oct 13 2017
 * @date last modification: Fri Jul 24 2020
 *
 * @brief  Facet synchronizer for parallel simulations with cohesive elments
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "element_synchronizer.hh"
#include "fe_engine.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_FACET_SYNCHRONIZER_HH_
#define AKANTU_FACET_SYNCHRONIZER_HH_

namespace akantu {

class FacetSynchronizer : public ElementSynchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  FacetSynchronizer(Mesh & mesh,
                    const ElementSynchronizer & element_synchronizer,
                    const ID & id = "facet_synchronizer");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// update distributed synchronizer after elements' insertion
  void
  updateDistributedSynchronizer(ElementSynchronizer & distributed_synchronizer,
                                DataAccessor<Element> & data_accessor,
                                const Mesh & mesh_cohesive);

protected:
  /// update elements list based on facets list
  void updateElementList(Array<Element> * elements,
                         const Array<Element> * facets,
                         const Mesh & mesh_cohesive);

  /// setup facet synchronization
  void
  setupFacetSynchronization(ElementSynchronizer & distributed_synchronizer);

  /// build send facet arrays
  void buildSendElementList(
      const Array<ElementTypeMapArray<UInt> *> & send_connectivity,
      const Array<ElementTypeMapArray<UInt> *> & recv_connectivity,
      const Array<ElementTypeMapArray<UInt> *> & temp_send_element);

  /// build recv facet arrays
  void buildRecvElementList(
      const Array<ElementTypeMapArray<UInt> *> & temp_recv_element);

  /// get facets' global connectivity for a list of elements
  template <GhostType ghost_facets>
  inline void getFacetGlobalConnectivity(
      const ElementSynchronizer & distributed_synchronizer,
      const ElementTypeMapArray<UInt> & rank_to_facet,
      const Array<Element> * elements,
      Array<ElementTypeMapArray<UInt> *> & connectivity,
      Array<ElementTypeMapArray<UInt> *> & facets);

  /// initialize ElementTypeMap containing correspondance between
  /// facets and processors
  void initRankToFacet(ElementTypeMapArray<UInt> & rank_to_facet);

  /// find which processor a facet is assigned to
  void buildRankToFacet(ElementTypeMapArray<UInt> & rank_to_facet,
                        const Array<Element> * elements);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  ElementTypeMapArray<UInt> facet_to_rank;
};

} // namespace akantu

#include "facet_synchronizer_inline_impl.hh"

#endif /* AKANTU_FACET_SYNCHRONIZER_HH_ */
