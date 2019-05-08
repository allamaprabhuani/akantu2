/**
 * @file contact_detection.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed Sep 12 2018
 * @date last modification: Tue Oct 23 2018
 *
 * @brief  Mother class for all detection algorithms
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_memory.hh"
#include "aka_grid_dynamic.hh"
#include "aka_bbox.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "fe_engine.hh"
#include "parsable.hh"
#include "element_group.hh"
#include "contact_element.hh"
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_DETECTOR_HH__
#define __AKANTU_CONTACT_DETECTOR_HH__


namespace akantu {

enum class Surface {
  master,
  slave
};
  
/* -------------------------------------------------------------------------- */ 

class ContactDetector :
    private Memory, public Parsable {

  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructors                                                  */
  /* ------------------------------------------------------------------------ */
public:
  
  ContactDetector(Mesh &, const ID & id = "contact_detector",
		  UInt memory_id = 0);

  ContactDetector(Mesh &, Array<Real> & positions,
		  const ID & id = "contact_detector",
		  UInt memory_id = 0);
        
  ~ContactDetector() = default;

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// performs the all search steps 
  void search(std::map<UInt, ContactElement> &);

  /// performs global spatial search to construct spatial grids
  void globalSearch(SpatialGrid<UInt> &, SpatialGrid<UInt> &);

  ///  performs local search to find closet master node to a slave node
  void localSearch(SpatialGrid<UInt> &, SpatialGrid<UInt> &);

  /// constructs contact map for a given pair of slave and master node
  void constructContactMap(std::map<UInt, ContactElement> &);
  
private:
  /// reads the input file to get contact detection options
  void parseSection();

  /// extracts vectors which forms the plane of element
  void vectorsAlongElement(const Element & /* element id     */,
			   Matrix<Real> &  /* vectors matrix */);

  /// computes orthogonal projection on master elements
  void computeOrthogonalProjection(const UInt &           /* slave node */,
				   const Array<Element> & /* master elements */,
				   Array<Real> &          /* normals */,
				   Array<Real> &          /* gaps */,
				   Array<Real> &          /* projections */);
 
  /// computes normal on an element
  void computeNormalOnElement(const Element & /* element id    */,
			      Vector<Real> &  /* normal vector */);

  /// computes tangents on a given natural coordinate
  void computeTangentsOnElement(const Element &, Vector<Real> &,
				Matrix<Real> &);

  /// computes projection of a query point on an element
  void computeProjectionOnElement(const Element &      /* element */,
				  const Vector<Real> & /* normal */,
				  const Vector<Real> & /* query */,
				  Vector<Real> &       /* projection
							  */,
				  Vector<Real> & /* real_projection */);

  /// computes natural projection of a real projection
  void computeNaturalProjection(const Element & /* element     */,
				Vector<Real> &  /* real projection  */,
				Vector<Real> &  /* natural projection */);

  /// compute normal projection of slave coord on a given element
  void normalProjection(const Element & el, const Vector<Real> & slave_coord,
			Vector<Real> & natural_coord, Real & tolerance);

  /* ------------------------------------------------------------------------ */
  /* Inline Methods                                                           */
  /* ------------------------------------------------------------------------ */
public:
  /// checks whether the natural projection is valid or not
  inline bool checkValidityOfProjection(Vector<Real> & );
  
  /// extracts the coordinates of an element
  inline void coordinatesOfElement(const Element & , Matrix<Real> & );

  /// computes the optimal cell size for grid 
  inline void computeCellSpacing(Vector<Real> & );

  /// constructs a grid containing nodes lying within bounding box
  inline void constructGrid(SpatialGrid<UInt> &, BBox &, const Array<UInt> &);

  /// constructs the bounding box based on nodes list
  inline void constructBoundingBox(BBox &, const Array<UInt> &);

  /// get the surface id
  template<Surface id>
  inline std::string getSurfaceId();
   
  /// set the surface id
  template<Surface id>
  inline void setSurfaceId(const std::string);
  
  /// computes the maximum in radius for a given mesh 
  inline void computeMaximalDetectionDistance();

  /// constructs the connectivity for a contact element
  inline Vector<UInt> constructConnectivity(UInt &, const Element &);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// returns the maximum detection distance
  AKANTU_GET_MACRO(MaximumDetectionDistance, max_dd, Real);

  /// sets the maximum detection distance
  AKANTU_SET_MACRO(MaximumDetectionDistance, max_dd, Real);

  /// returns the bounding box extension
  AKANTU_GET_MACRO(MaximumBoundingBox, max_bb, Real);

  /// sets the bounding box extension
  AKANTU_SET_MACRO(MaximumBoundingBox, max_bb, Real);

  /// get the contact pairs
  AKANTU_GET_MACRO(ContactPairs, contact_pairs, Array<UInt>);
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// maximal detection distance for grid spacing
  Real max_dd;

  /// maximal bounding box extension
  Real max_bb;
  
  /// Mesh
  Mesh & mesh;

  /// dimension of the model
  UInt spatial_dimension{0};
  
  /// map to contain ids for surfaces
  std::map<Surface, std::string> surfaces; 

  /// contact pair slave node to closet master node
  Array<UInt> contact_pairs;
  
  /// contains the updated positions of the nodes
  Array<Real> & positions;

  /// type of detection explicit/implicit 
  DetectionType detection_type;
};
  
} // namespace akantu

#include "contact_detector_inline_impl.cc"

#endif /* __AKANTU_CONTACT_DETECTOR_HH__ */
