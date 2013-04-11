/**
 * @file   facet_synchronizer_inline_impl.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Thu Mar 28 13:08:23 2013
 *
 * @brief  facet synchronizer inline implementation
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
template<GhostType ghost_facets>
inline void FacetSynchronizer::getFacetBarycentersPerElement(const ByElementTypeUInt & rank_to_facet,
							     const Array<Element> * elements,
							     std::list<ElementBarycenter> * facet_elbary) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh_facets.getSpatialDimension();
  ElementBarycenter current_elbary(spatial_dimension);

  /// loop on every processor
  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == rank) continue;
    std::list<ElementBarycenter> & elbary = facet_elbary[p];
    const Array<Element> & elem = elements[p];

    UInt nb_element = elem.getSize();

    /// loop on every send/recv element
    for (UInt el = 0; el < nb_element; ++el) {
      ElementType type = elem(el).type;
      GhostType gt = elem(el).ghost_type;
      UInt el_index = elem(el).element;

      const Array<Element> & facet_to_element =
	mesh_facets.getSubelementToElement(type, gt);
      UInt nb_facets_per_element = Mesh::getNbFacetsPerElement(type);
      ElementType facet_type = Mesh::getFacetType(type);

      /// loop on every facet of the element
      for (UInt f = 0; f < nb_facets_per_element; ++f) {

	current_elbary.elem = facet_to_element(el_index, f);
	UInt facet_index = current_elbary.elem.element;
	GhostType facet_gt = current_elbary.elem.ghost_type;

	const Array<UInt> & t_to_f = rank_to_facet(facet_type, facet_gt);

	/// exclude not ghost facets, facets assigned to other
	/// processors
	if (facet_gt != ghost_facets) continue;
	if ((facet_gt == _ghost) && (t_to_f(facet_index) != p)) continue;

	/// exclude facets on the boundary
	if (ghost_facets == _ghost) {
	  const Array<std::vector<Element> > & element_to_facet =
	    mesh_facets.getElementToSubelement(facet_type, facet_gt);
	  if (element_to_facet(facet_index)[1] == ElementNull) continue;
	}

	/// store barycenter's coordinates
	mesh_facets.getBarycenter(facet_index,
				  facet_type,
				  current_elbary.bary.storage(),
				  facet_gt);
	elbary.push_back(current_elbary);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}
