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
  ///
  void search(std::map<UInt, ContactElement> &);

  /// computes orthogonal projection on master elements
  void computeOrthogonalProjection(const UInt &           /* slave node */,
				   const Array<Element> & /* master elements */,
				   Array<Real> &          /* normals */,
				   Array<Real> &          /* projections */);

private:
  /// reads the input file to get contact detection options
  void parseSection();
     
  /// performs global spatial search
  void globalSearch(std::map<UInt, ContactElement> &);

  ///  performs local search to create contact element
  /// TODO: templated function typename
  void localSearch(SpatialGrid<UInt> &, SpatialGrid<UInt> &,
		   std::map<UInt, ContactElement> &);

  /// constructs a grid containing nodes lying within bounding box
  /// TODO : templated fucntion to created template Spatial Grid
  void constructGrid(SpatialGrid<UInt> &, BBox &, const Array<UInt> &);

  /// constructs the bounding box based on nodes list
  void constructBoundingBox(BBox &, const Array<UInt> &);

  /// computes the optimal cell size for grid 
  void computeCellSpacing(Vector<Real> &);

  /// computes the maximum in radius for a given mesh 
  void getMaximalDetectionDistance();

  /// extracts the coordinates of an element
  void coordinatesOfElement(const Element & /* element id  */,
			    Matrix<Real> &  /* coordinates */);

  /// extracts vectors which forms the plane of element
  void vectorsAlongElement(const Element & /* element id     */,
			   Matrix<Real> &  /* vectors matrix */);
  
  /// computes normal on an element
  void computeNormalOnElement(const Element & /* element id    */,
			      Vector<Real> &  /* normal vector */);

  /// computes projection of a query point on an element
  void computeProjectionOnElement(const Element &      /* element */,
				  const Vector<Real> & /* normal */,
				  const Vector<Real> & /* query */,
				  Vector<Real> &       /* projection
							  */);
  /// checks for the validity of a projection
  bool isValidProjection(const Element & /* element     */,
			 Vector<Real> &  /* projection  */);
  
  
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
  
  /// contains the updated positions of the nodes
  Array<Real> & positions;

  /// type of detection extrinisic/intrinsic 
  ContactDetectorType detection_type;
  
};
  
} // namespace akantu


#endif /* __AKANTU_CONTACT_DETECTOR_HH__ */
