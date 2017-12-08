/**
 * @file   cohesive_element_inserter.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Dec 04 2013
 * @date last modification: Fri Oct 02 2015
 *
 * @brief  Cohesive element inserter
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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
#include "data_accessor.hh"
#include "mesh_utils.hh"
#include "parsable.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COHESIVE_ELEMENT_INSERTER_HH__
#define __AKANTU_COHESIVE_ELEMENT_INSERTER_HH__

namespace akantu {
class GlobalIdsUpdater;
class FacetSynchronizer;
} // akantu

namespace akantu {

class CohesiveElementInserter : public DataAccessor<Element>, public Parsable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  CohesiveElementInserter(Mesh & mesh, bool is_extrinsic = false,
                          const ID & id = "cohesive_element_inserter");

  ~CohesiveElementInserter() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// init function
  void init(bool is_extrinsic);

  /// set range limitation for intrinsic cohesive element insertion
  void setLimit(SpacialDirection axis, Real first_limit, Real second_limit);

  /// insert intrinsic cohesive elements in a predefined range
  UInt insertIntrinsicElements();

  /// preset insertion of intrinsic cohesive elements along
  /// a predefined group of facet and assign them a defined material index.
  /// insertElement() method has to be called to finalize insertion.
  void insertIntrinsicElements(const std::string & physname,
                               UInt material_index);

  /// insert extrinsic cohesive elements (returns the number of new
  /// cohesive elements)
  UInt insertElements(bool only_double_facets = false);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// limit check facets to match given insertion limits
  void limitCheckFacets();

  /// init parallel variables
  void initParallel(ElementSynchronizer & element_synchronizer);

protected:
  /// init facets check
  void initFacetsCheck();

  /// update facet insertion arrays after facets doubling
  void updateInsertionFacets();

  /// update nodes type and global ids for parallel simulations
  /// (locally, within each processor)
  UInt updateGlobalIDs(NewNodesEvent & node_event);

  /// synchronize the global ids among the processors in parallel simulations
  void synchronizeGlobalIDs(NewNodesEvent & node_event);

  /// update nodes type
  void updateNodesType(Mesh & mesh, NewNodesEvent & node_event);

  /// functions for parallel communications
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  template <bool pack_mode>
  inline void
  packUnpackGlobalConnectivity(CommunicationBuffer & buffer,
                               const Array<Element> & elements) const;

  template <bool pack_mode>
  inline void
  packUnpackGroupedInsertionData(CommunicationBuffer & buffer,
                                 const Array<Element> & elements) const;
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO_NOT_CONST(InsertionFacetsByElement, insertion_facets,
                             ElementTypeMapArray<bool> &);

  AKANTU_GET_MACRO(InsertionFacetsByElement, insertion_facets,
                   const ElementTypeMapArray<bool> &);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(InsertionFacets, insertion_facets, bool);

  AKANTU_GET_MACRO(CheckFacets, check_facets, const ElementTypeMapArray<bool> &);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(CheckFacets, check_facets, bool);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(CheckFacets, check_facets, bool);
  AKANTU_GET_MACRO(MeshFacets, mesh_facets, const Mesh &);
  AKANTU_GET_MACRO_NOT_CONST(MeshFacets, mesh_facets, Mesh &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// object id
  ID id;

  /// main mesh where to insert cohesive elements
  Mesh & mesh;

  /// mesh containing facets
  Mesh & mesh_facets;

  /// list of facets where to insert elements
  ElementTypeMapArray<bool> insertion_facets;

  /// limits for element insertion
  Matrix<Real> insertion_limits;

  /// list of groups to consider for insertion, ignored if empty
  std::vector<ID> physical_groups;

  /// vector containing facets in which extrinsic cohesive elements can be
  /// inserted
  ElementTypeMapArray<bool> check_facets;

  /// global connectivity ids updater
  std::unique_ptr<GlobalIdsUpdater> global_ids_updater;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const CohesiveElementInserter & _this) {
  _this.printself(stream);
  return stream;
}

class CohesiveNewNodesEvent : public NewNodesEvent {
public:
  CohesiveNewNodesEvent() = default;
  ~CohesiveNewNodesEvent() override = default;

  AKANTU_GET_MACRO_NOT_CONST(OldNodesList, old_nodes, Array<UInt> &);
  AKANTU_GET_MACRO(OldNodesList, old_nodes, const Array<UInt> &);

private:
  Array<UInt> old_nodes;
};

} // akantu

#include "cohesive_element_inserter_inline_impl.cc"

#endif /* __AKANTU_COHESIVE_ELEMENT_INSERTER_HH__ */
