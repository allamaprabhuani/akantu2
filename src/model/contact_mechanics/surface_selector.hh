/**
 * @file   surface_selector.hh
 *
 * @author Mohit Pundir <pundir.pundir@epfl.ch>
 *
 * @date creation: Fri Jun 21 2019
 * @date last modification: Fri Jun 21 2019
 *
 * @brief  Node selectors for contact detection
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh.hh"
#include "mesh_utils.hh"
#include "parsable.hh"
#if defined(AKANTU_COHESIVE_ELEMENT)
#include "cohesive_element_inserter.hh"
#endif

/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SURFACE_SELECTOR_HH__
#define __AKANTU_SURFACE_SELECTOR_HH__

namespace akantu {
class Model;
class GlobalIdsUpdater;
} // namespace akantu
namespace akantu {

/**
 * main class to assign surfaces for contact detection
 */
class SurfaceSelector : public MeshEventHandler, public Parsable {
  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructor                                                   */
  /* ------------------------------------------------------------------------ */
public:
  SurfaceSelector(Mesh & mesh);

  virtual ~SurfaceSelector() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual Array<UInt> & getMasterList() = 0;

  virtual Array<UInt> & getSlaveList() = 0;

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  Mesh & mesh;
};

/* -------------------------------------------------------------------------- */

class PhysicalSurfaceSelector : public SurfaceSelector {
  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructor                                                   */
  /* ------------------------------------------------------------------------ */
public:
  PhysicalSurfaceSelector(Mesh & mesh);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  Array<UInt> & getMasterList() override;

  Array<UInt> & getSlaveList() override;

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  std::string master;
  std::string slave;
};

/* -------------------------------------------------------------------------- */
#if defined(AKANTU_COHESIVE_ELEMENT)
class CohesiveSurfaceSelector : public SurfaceSelector {
  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructor                                                   */
  /* ------------------------------------------------------------------------ */
public:
  CohesiveSurfaceSelector(Mesh & mesh);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  void onNodesAdded(const Array<UInt> & nodes_list,
                    const NewNodesEvent & event) override;

  void filterBoundaryElements(Array<Element> & elements,
                              Array<Element> & boundary_elements);

public:
  Array<UInt> & getMasterList() override;

  Array<UInt> & getSlaveList() override;

  /* ------------------------------------------------------------------------
   */
  /*                                                                          */
  /* ------------------------------------------------------------------------
   */
  AKANTU_GET_MACRO_NOT_CONST(NewNodesList, new_nodes_list, Array<UInt> &);
  AKANTU_GET_MACRO(NewNodesList, new_nodes_list, const Array<UInt> &);

  /* ------------------------------------------------------------------------
   */
  /* Members */
  /* ------------------------------------------------------------------------
   */
protected:
  Mesh & mesh_facets;

  Array<UInt> new_nodes_list;
};

class AllSurfaceSelector : public SurfaceSelector {
  /* ------------------------------------------------------------------------
   */
  /* Constructor/Destructor */
  /* ------------------------------------------------------------------------
   */
public:
  AllSurfaceSelector(Mesh & mesh);

  /* ------------------------------------------------------------------------
   */
  /* Methods */
  /* ------------------------------------------------------------------------
   */
protected:
  void onNodesAdded(const Array<UInt> & nodes_list,
                    const NewNodesEvent & event) override;

  void filterBoundaryElements(Array<Element> & elements,
                              Array<Element> & boundary_elements);

public:
  Array<UInt> & getMasterList() override;

  Array<UInt> & getSlaveList() override;

  /* ------------------------------------------------------------------------
   */
  /*                                                                          */
  /* ------------------------------------------------------------------------
   */
  AKANTU_GET_MACRO_NOT_CONST(NewNodesList, new_nodes_list, Array<UInt> &);
  AKANTU_GET_MACRO(NewNodesList, new_nodes_list, const Array<UInt> &);

  /* ------------------------------------------------------------------------
   */
  /* Members */
  /* ------------------------------------------------------------------------
   */
protected:
  std::string master;

  std::string slave;

  Mesh & mesh_facets;

  Array<UInt> new_nodes_list;
};

#endif

} // namespace akantu

#endif /*  __AKANTU_SURFACE_SELECTOR_HH__  */
