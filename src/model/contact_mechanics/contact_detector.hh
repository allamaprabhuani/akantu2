/**
 * @file   contact_detector.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Thu Jun 24 2021
 *
 * @brief  Mother class for all detection algorithms
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
#include "aka_bbox.hh"
#include "aka_common.hh"
#include "aka_grid_dynamic.hh"
#include "contact_element.hh"
#include "element_class.hh"
#include "element_group.hh"
#include "fe_engine.hh"
#include "geometry_utils.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "parsable.hh"
#include "surface_selector.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_DETECTOR_HH__
#define __AKANTU_CONTACT_DETECTOR_HH__

namespace akantu {

enum class Surface { master, slave };

/* -------------------------------------------------------------------------- */

class ContactDetector : public Parsable {

  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructors                                                  */
  /* ------------------------------------------------------------------------ */
public:
  ContactDetector(Mesh & /*mesh*/, const ID & id = "contact_detector");

  ContactDetector(Mesh & /*mesh*/, Array<Real> positions,
                  const ID & id = "contact_detector");

  ~ContactDetector() override = default;

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// performs all search steps
  void search(Array<ContactElement> & elements, Array<Real> & gaps,
              Array<Real> & normals, Array<Real> & tangents,
              Array<Real> & projections);

  /// performs global spatial search to construct spatial grids
  void globalSearch(SpatialGrid<UInt> & /*slave_grid*/,
                    SpatialGrid<UInt> & /*master_grid*/);

  ///  performs local search to find closet master node to a slave node
  void localSearch(SpatialGrid<UInt> & /*slave_grid*/,
                   SpatialGrid<UInt> & /*master_grid*/);

  /// create contact elements
  void createContactElements(Array<ContactElement> & elements,
                             Array<Real> & gaps, Array<Real> & normals,
                             Array<Real> & tangents, Array<Real> & projections);

private:
  /// reads the input file to get contact detection options
  void parseSection(const ParserSection & section) override;

  /* ------------------------------------------------------------------------ */
  /* Inline Methods                                                           */
  /* ------------------------------------------------------------------------ */
public:
  /// checks whether the natural projection is valid or not
  inline bool checkValidityOfProjection(Vector<Real> & projection) const;

  /// extracts the coordinates of an element
  inline void coordinatesOfElement(const Element & el,
                                   Matrix<Real> & coords) const;

  /// computes the optimal cell size for grid
  inline void computeCellSpacing(Vector<Real> & spacing) const;

  /// constructs a grid containing nodes lying within bounding box
  inline void constructGrid(SpatialGrid<UInt> & grid, BBox & bbox,
                            const Array<UInt> & nodes_list) const;

  /// constructs the bounding box based on nodes list
  inline void constructBoundingBox(BBox & bbox,
                                   const Array<UInt> & nodes_list) const;

  /// computes the maximum in radius for a given mesh
  inline void computeMaximalDetectionDistance();

  /// constructs the connectivity for a contact element
  inline Vector<UInt> constructConnectivity(UInt & slave,
                                            const Element & master) const;

  /// computes normal on an element
  inline void computeNormalOnElement(const Element & element,
                                     Vector<Real> & normal) const;

  /// extracts vectors which forms the plane of element
  inline void vectorsAlongElement(const Element & el,
                                  Matrix<Real> & vectors) const;

  /// computes the gap between slave and its projection on master
  /// surface
  inline Real computeGap(const Vector<Real> & slave,
                         const Vector<Real> & master) const;

  /// filter boundary elements
  inline void filterBoundaryElements(const Array<Element> & elements,
                                     Array<Element> & boundary_elements) const;

  /// checks whether self contact condition leads to a master element
  /// which is closet but not orthogonally opposite to slave surface
  inline bool isValidSelfContact(const UInt & slave_node, const Real & gap,
                                 const Vector<Real> & normal) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the mesh
  AKANTU_GET_MACRO(Mesh, mesh, Mesh &)

  /// returns the maximum detection distance
  AKANTU_GET_MACRO(MaximumDetectionDistance, max_dd, Real);
  AKANTU_SET_MACRO(MaximumDetectionDistance, max_dd, Real);

  /// returns the bounding box extension
  AKANTU_GET_MACRO(MaximumBoundingBox, max_bb, Real);
  AKANTU_SET_MACRO(MaximumBoundingBox, max_bb, Real);

  /// returns the minimum detection distance
  AKANTU_GET_MACRO(MinimumDetectionDistance, min_dd, Real);
  AKANTU_SET_MACRO(MinimumDetectionDistance, min_dd, Real);

  AKANTU_GET_MACRO_NOT_CONST(Positions, positions, Array<Real> &);
  AKANTU_SET_MACRO(Positions, positions, Array<Real>);

  AKANTU_GET_MACRO_NOT_CONST(SurfaceSelector, *surface_selector,
                             SurfaceSelector &);
  AKANTU_SET_MACRO(SurfaceSelector, surface_selector,
                   std::shared_ptr<SurfaceSelector>);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// maximal detection distance for grid spacing
  Real max_dd;

  /// minimal detection distance for grid spacing
  Real min_dd;

  /// maximal bounding box extension
  Real max_bb;

  /// tolerance for finding natural projection
  Real projection_tolerance;

  /// iterations for finding natural projection
  UInt max_iterations;

  /// tolerance for extending a master elements on all sides
  Real extension_tolerance;

  /// Mesh
  Mesh & mesh;

  /// dimension of the model
  UInt spatial_dimension{0};

  /// node selector for selecting master and slave nodes
  std::shared_ptr<SurfaceSelector> surface_selector;

  /// contact pair slave node to closet master node
  std::vector<std::pair<UInt, UInt>> contact_pairs;

  /// contains the updated positions of the nodes
  Array<Real> positions;

  /// type of detection explicit/implicit
  DetectionType detection_type;
};

} // namespace akantu

#include "contact_detector_inline_impl.cc"

#endif /* __AKANTU_CONTACT_DETECTOR_HH__ */
