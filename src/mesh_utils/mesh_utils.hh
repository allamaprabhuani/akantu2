/**
 * @file   mesh_utils.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Fri Aug 20 12:19:44 2010
 *
 * @brief  All mesh utils necessary for various tasks
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

#include "aka_common.hh"
#include "mesh.hh"
#include "aka_csr.hh"
/* -------------------------------------------------------------------------- */

#include <vector>

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MESH_UTILS_HH__
#define __AKANTU_MESH_UTILS_HH__

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

class MeshUtils {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MeshUtils();
  virtual ~MeshUtils();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// build map from nodes to elements
  static void buildNode2Elements(const Mesh & mesh, 
				 CSR<UInt> & node_to_elem,
				 UInt spatial_dimension = _all_dimensions);
  static void buildNode2Elements(const Mesh & mesh,
				 CSR<Element> & node_to_elem,
				 UInt spatial_dimension = _all_dimensions);

  /// build map from nodes to elements for a specific element type
  static void buildNode2ElementsByElementType(const Mesh & mesh,
					      CSR<UInt> & node_to_elem,
					      const ElementType & type,
					      const GhostType & ghost_type = _not_ghost);

  /// build facets elements on boundary
  static void buildFacets(Mesh & mesh);

  /// build facets elements : boundary and internals (for serial simulations)
  static void buildAllFacets(Mesh & mesh,
			     Mesh & mesh_facets);

  /// build facets elements : boundary and internals (for parallel simulations)
  static void buildAllFacetsParallel(Mesh & mesh,
				     Mesh & mesh_facets,
				     ByElementTypeUInt & prank_to_element);

  /// build facets for a given spatial dimension
  static void buildFacetsDimension(Mesh & mesh,
				   Mesh & mesh_facets,
				   bool boundary_only,
				   UInt dimension,
				   ByElementTypeUInt & prank_to_element);

  /// build normal to some elements
  //  static void buildNormals(Mesh & mesh, UInt spatial_dimension=0);

  /// take  the local_connectivity  array  as  the array  of  local and  ghost
  /// connectivity, renumber the nodes and set the connectivity of the mesh
  static void renumberMeshNodes(Mesh & mesh,
				UInt * local_connectivities,
				UInt nb_local_element,
				UInt nb_ghost_element,
				ElementType type,
				Array<UInt> & old_nodes);

//  static void setUIntData(Mesh & mesh, UInt * data, UInt nb_tags, const ElementType & type);

  /// Detect closed surfaces of the mesh and save the surface id
  /// of the surface elements in the array surface_id
  static void buildSurfaceID(Mesh & mesh);

  /// compute pbc pair for on given direction
  static void computePBCMap(const Mesh & mymesh,const UInt dir,
		     std::map<UInt,UInt> & pbc_pair);
  /// compute pbc pair for a surface pair
  static void computePBCMap(const Mesh & mymesh,
			    const std::pair<Surface, Surface> & surface_pair,
			    const ElementType type,
			    std::map<UInt,UInt> & pbc_pair);


  // /// tweak mesh connectivity to activate pbc
  // static void tweakConnectivityForPBC(Mesh & mesh,
  // 				      bool flag_x,
  // 				      bool flag_y = false,
  // 				      bool flag_z = false);

  /// create a multimap of nodes per surfaces
  static void buildNodesPerSurface(const Mesh & mesh, CSR<UInt> & nodes_per_surface);

  /// function to print the contain of the class
  //  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// remove not connected nodes /!\ this functions renumbers the nodes.
  static void purifyMesh(Mesh & mesh);

  /// function to insert intrinsic cohesive elements in a zone
  /// delimited by provided limits
  static void insertIntrinsicCohesiveElementsInArea(Mesh & mesh,
						    const Array<Real> & limits);

  /// function to insert intrinsic cohesive elements on the selected
  /// facets
  static void insertIntrinsicCohesiveElements(Mesh & mesh,
					      Mesh & mesh_facets,
					      ElementType type_facet,
					      const Array<bool> & facet_insertion);

  /// function to insert cohesive elements on the selected facets
  static void insertCohesiveElements(Mesh & mesh,
				     Mesh & mesh_facets,
				     ByElementTypeArray<bool> & facet_insertion,
				     const bool extrinsic);

  /// fill the subelement to element and the elements to subelements data
  static void fillElementToSubElementsData(Mesh & mesh);

  /// flip facets based on global connectivity
  static void flipFacets(Mesh & mesh_facets,
			 const ByElementTypeUInt & global_connectivity,
			 GhostType gt_facet);

  /// provide list of elements around a node and check if a given
  /// facet is reached
  template <bool third_dim_points>
  static bool findElementsAroundSubfacet(const Mesh & mesh,
					 const Mesh & mesh_facets,
					 const Element & starting_element,
					 const Element & end_facet,
					 const Vector<UInt> & subfacet_connectivity,
					 std::vector<Element> & elem_list,
					 std::vector<Element> & facet_list,
					 std::vector<Element> * subfacet_list = NULL);

  /// function to check if a node belongs to a given element
  static inline bool hasElement(const Array<UInt> & connectivity,
				const Element & el,
				const Vector<UInt> & nodes);

  /// reset facet to double arrays
  static void resetFacetToDouble(Mesh & mesh_facets);

private:

  /// match pairs that are on the associated pbc's
  static void matchPBCPairs(const Mesh & mymesh,
			    const UInt dir,
			    std::vector<UInt> & selected_left,
			    std::vector<UInt> & selected_right,
			    std::map<UInt,UInt> & pbc_pair);

  /// function used by all the renumbering functions
  static void renumberNodesInConnectivity(UInt * list_nodes,
					  UInt nb_nodes,
					  std::map<UInt, UInt> & renumbering_map);

  /// update facet_to_subfacet
  static void updateFacetToSubfacet(Mesh & mesh_facets,
				    ElementType type_subfacet,
				    GhostType gt_subfacet,
				    bool facet_mode);

  /// update subfacet_to_facet
  static void updateSubfacetToFacet(Mesh & mesh_facets,
				    ElementType type_subfacet,
				    GhostType gt_subfacet,
				    bool facet_mode);

  /// function to double a given facet and update the list of doubled
  /// nodes
  static void doubleFacet(Mesh & mesh,
			  Mesh & mesh_facets,
			  UInt facet_dimension,
			  Array<UInt> & doubled_nodes,
			  bool facet_mode);

  /// function to double a subfacet given start and end index for
  /// local facet_to_subfacet vector, and update the list of doubled
  /// nodes
  template <UInt spatial_dimension>
  static void doubleSubfacet(Mesh & mesh,
			     Mesh & mesh_facets,
			     Array<UInt> & doubled_nodes);

  /// double a node
  static void doubleNodes(Mesh & mesh,
			  const Array<UInt> & old_nodes,
			  Array<UInt> & doubled_nodes);

  /// fill facet_to_double array in the mesh
  static bool updateFacetToDouble(Mesh & mesh_facets,
				  const ByElementTypeArray<bool> & facet_insertion);

  /// find subfacets to be doubled
  template <bool subsubfacet_mode>
  static void findSubfacetToDouble(Mesh & mesh, Mesh & mesh_facets);

  /// double facets (points) in 1D
  static void doublePointFacet(Mesh & mesh,
			       Mesh & mesh_facets,
			       Array<UInt> & doubled_nodes);

  /// update cohesive element data
  static void updateCohesiveData(Mesh & mesh,
				 Mesh & mesh_facets,
				 NewElementsEvent & element_event);

  /// update elemental connectivity after doubling a node
  inline static void updateElementalConnectivity(Mesh & mesh,
						 UInt old_node,
						 UInt new_node,
						 const std::vector<Element> & element_list,
						 const std::vector<Element> * facet_list = NULL);

  /// double middle nodes if facets are _segment_3
  template <bool third_dim_segments>
  static void updateQuadraticSegments(Mesh & mesh,
				      Mesh & mesh_facets,
				      ElementType type_facet,
				      GhostType gt_facet,
				      Array<UInt> & doubled_nodes);

  /// update nodes type and global ids for parallel simulations
  static UInt updateGlobalIDs(Mesh & mesh,
			      Mesh & mesh_facets,
			      const Array<UInt> & doubled_nodes);

  /// remove elements on a vector
  inline static bool removeElementsInVector(const std::vector<Element> & elem_to_remove,
					    std::vector<Element> & elem_list);

  static void sortElements(std::vector<Element> & elements, const Vector<UInt> facet,
                           const Mesh & mesh, const Mesh & mesh_facets,
                           const ByElementTypeReal & barycenters);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

};


class ElementSorter {
public:
  ElementSorter(const ElementSorter & e) : atan2(e.atan2) {}

  ElementSorter(std::map<Element, Real, CompElementLess> & atan2) : atan2(atan2) {}

  inline bool operator()(const Element & first, const Element & second);
private:
  std::map<Element, Real, CompElementLess> & atan2;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "mesh_utils_inline_impl.cc"
#endif

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const MeshUtils & _this)
// {
//   _this.printself(stream);
//   return stream;
// }

__END_AKANTU__
#endif /* __AKANTU_MESH_UTILS_HH__ */
