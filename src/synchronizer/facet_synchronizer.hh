/**
 * Copyright (©) 2017-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
      const Array<ElementTypeMapArray<Idx> *> & send_connectivity,
      const Array<ElementTypeMapArray<Idx> *> & recv_connectivity,
      const Array<ElementTypeMapArray<Idx> *> & temp_send_element);

  /// build recv facet arrays
  void buildRecvElementList(
      const Array<ElementTypeMapArray<Idx> *> & temp_recv_element);

  /// get facets' global connectivity for a list of elements
  template <GhostType ghost_facets>
  inline void getFacetGlobalConnectivity(
      const ElementSynchronizer & distributed_synchronizer,
      const ElementTypeMapArray<Idx> & rank_to_facet,
      const Array<Element> * elements,
      Array<ElementTypeMapArray<Idx> *> & connectivity,
      Array<ElementTypeMapArray<Idx> *> & facets);

  /// initialize ElementTypeMap containing correspondance between
  /// facets and processors
  void initRankToFacet(ElementTypeMapArray<Idx> & rank_to_facet);

  /// find which processor a facet is assigned to
  void buildRankToFacet(ElementTypeMapArray<Idx> & rank_to_facet,
                        const Array<Element> * elements);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  ElementTypeMapArray<Int> facet_to_rank;
};

} // namespace akantu

#include "facet_synchronizer_inline_impl.hh"

#endif /* AKANTU_FACET_SYNCHRONIZER_HH_ */
