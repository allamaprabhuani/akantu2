/**
 * Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh.hh"
#include "mesh_utils.hh"
#include "parsable.hh"
#if defined(AKANTU_COHESIVE_ELEMENT)
#include "cohesive_element_inserter.hh"
#endif

/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SURFACE_SELECTOR_HH_
#define AKANTU_SURFACE_SELECTOR_HH_

namespace akantu {
class Model;
class GlobalIdsUpdater;
} // namespace akantu

namespace akantu {

/**
 * main class to assign surfaces for contact detection
 */
class SurfaceSelector : public MeshEventHandler, public Parsable {
public:
  SurfaceSelector(Mesh & mesh);
  ~SurfaceSelector() override = default;

public:
  virtual Array<Idx> & getMasterList() { AKANTU_TO_IMPLEMENT(); }
  virtual Array<Idx> & getSlaveList() { AKANTU_TO_IMPLEMENT(); }

protected:
  Mesh & mesh;
};

/* -------------------------------------------------------------------------- */
/**
 * class that selects contact surface from physical names
 */
class PhysicalSurfaceSelector : public SurfaceSelector {
public:
  PhysicalSurfaceSelector(Mesh & mesh);

public:
  Array<Idx> & getMasterList() override;
  Array<Idx> & getSlaveList() override;

protected:
  std::string master;
  std::string slave;
};

/* -------------------------------------------------------------------------- */
/**
 * class that selects contact surface from cohesive elements
 */
#if defined(AKANTU_COHESIVE_ELEMENT)
class CohesiveSurfaceSelector : public SurfaceSelector {
public:
  CohesiveSurfaceSelector(Mesh & mesh);

protected:
  void onElementsAdded(const Array<Element> & element_list,
                       const NewElementsEvent & event) override;

  void onNodesAdded(const Array<Idx> & nodes_list,
                    const NewNodesEvent & event) override;

public:
  Array<Idx> & getMasterList() override;
  Array<Idx> & getSlaveList() override;

  AKANTU_GET_MACRO_NOT_CONST(NewNodesList, new_nodes_list, Array<Idx> &);
  AKANTU_GET_MACRO(NewNodesList, new_nodes_list, const Array<Idx> &);

protected:
  Mesh & mesh_facets;
  Array<Idx> new_nodes_list;
};

/* -------------------------------------------------------------------------- */
/**
 * class that selects contact surface from both cohesive elements and
 * physical names
 */
class AllSurfaceSelector : public SurfaceSelector {
public:
  AllSurfaceSelector(Mesh & mesh);

protected:
  void onElementsAdded(const Array<Element> & element_list,
                       const NewElementsEvent & event) override;

  void onNodesAdded(const Array<Int> & nodes_list,
                    const NewNodesEvent & event) override;

public:
  Array<Idx> & getMasterList() override;

  Array<Idx> & getSlaveList() override;

  AKANTU_GET_MACRO_NOT_CONST(NewNodesList, new_nodes_list, Array<Idx> &);
  AKANTU_GET_MACRO(NewNodesList, new_nodes_list, const Array<Idx> &);

protected:
  std::string master;
  std::string slave;
  Mesh & mesh_facets;
  Array<Idx> new_nodes_list;
};

#endif

} // namespace akantu

#endif /*  AKANTU_SURFACE_SELECTOR_HH_  */
