/**
 * @file   facet_synchronizer.hh
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Tue Mar 26 09:49:44 2013
 *
 * @brief  Facet synchronizer for parallel simulations with cohesive elments
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

#ifndef __AKANTU_FACET_SYNCHRONIZER_HH__
#define __AKANTU_FACET_SYNCHRONIZER_HH__

#include "distributed_synchronizer.hh"
#include "fem.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class FacetSynchronizer : public DistributedSynchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
protected:

  FacetSynchronizer(DistributedSynchronizer & distributed_synchronizer,
		    Mesh & mesh_facets,
		    SynchronizerID id = "facet_synchronizer",
		    MemoryID memory_id = 0);

// public:

//   virtual ~FacetSynchronizer() {};

  /* ------------------------------------------------------------------------ */
  /* Subclasses                                                               */
  /* ------------------------------------------------------------------------ */

protected:

  class ElementBarycenter {
  public:
    /// constructor/deconstructor
    ElementBarycenter(UInt spatial_dimension) : bary(spatial_dimension) {};
    ~ElementBarycenter() {};

    /// members
    Element elem;
    Vector<Real> bary;

    /// operator < for sorting
    inline bool operator<(const ElementBarycenter & x) {
      UInt size = this->bary.size();
      return std::lexicographical_compare(this->bary.storage(),
					  this->bary.storage() + size,
					  x.bary.storage(),
					  x.bary.storage() + size);
    };
  };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// get a distributed synchronizer and a facet mesh, and create the
  /// associated FacetSynchronizer
  static FacetSynchronizer *
  createFacetSynchronizer(DistributedSynchronizer & distributed_synchronizer,
			  Mesh & mesh_facets,
			  SynchronizerID id = "facet_synchronizer",
			  MemoryID memory_id = 0);

  /// update distributed synchronizer after elements' insertion
  void updateDistributedSynchronizer(DataAccessor & data_accessor,
				     const ByElementTypeUInt & cohesive_el_to_facet);

protected:

  /// update elements list based on facets list
  void updateElementList(Array<Element> * elements,
			 const Array<Element> * facets,
			 const ByElementTypeUInt & cohesive_el_to_facet);

  /// setup facet synchronization
  void setupFacetSynchronization();

  /// build send/recv facet arrays
  void buildFacetList(std::list<ElementBarycenter> * element_barycenter,
		      Array<Element> * final_facets);

  /// get facets' barycenter for a list of elements
  template<GhostType ghost_facets>
  inline void getFacetBarycentersPerElement(const ByElementTypeUInt & rank_to_facet,
					    const Array<Element> * elements,
					    std::list<ElementBarycenter> * facet_elbary);

  /// initialize ByElementType containing correspondance between
  /// facets and processors
  void initRankToFacet(ByElementTypeUInt & rank_to_facet);

  /// find which processor a facet is assigned to
  void buildRankToFacet(ByElementTypeUInt & rank_to_facet,
			const Array<Element> * elements);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  DistributedSynchronizer & distributed_synchronizer;

  Mesh & mesh_facets;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "facet_synchronizer_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_FACET_SYNCHRONIZER_HH__ */
