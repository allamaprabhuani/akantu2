/**
 * @file   material_selector.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Aug 13 10:51:44 2013
 *
 * @brief  class describing how to choose a material for a given element
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



#ifndef __AKANTU_MATERIAL_SELECTOR_HH__
#define __AKANTU_MATERIAL_SELECTOR_HH__

__BEGIN_AKANTU__

class SolidMechanicsModel;

/* -------------------------------------------------------------------------- */
class MaterialSelector {
public:
  MaterialSelector() : fallback_value(0) {}
  virtual ~MaterialSelector() {}
  virtual UInt operator()(const Element & element) {
    return fallback_value;
  }

  void setFallback(UInt f) { fallback_value = f; }
private:
  UInt fallback_value;
};

/* -------------------------------------------------------------------------- */
class DefaultMaterialSelector : public MaterialSelector {
public:
  DefaultMaterialSelector(const ByElementTypeUInt & element_index_by_material) :
    element_index_by_material(element_index_by_material) { }

  UInt operator()(const Element & element) {
    try {
      DebugLevel dbl = debug::getDebugLevel();
      debug::setDebugLevel(dblError);
      UInt mat = element_index_by_material(element.type, element.ghost_type)(element.element, 0);
      debug::setDebugLevel(dbl);
      return mat;
    } catch (...) {
      return MaterialSelector::operator()(element);
    }
  }

private:
  const ByElementTypeUInt & element_index_by_material;
};

/* -------------------------------------------------------------------------- */
template<typename T>
class MeshDataMaterialSelector : public MaterialSelector {
public:
  MeshDataMaterialSelector(const std::string & name, const SolidMechanicsModel & model) : mesh_data(name), model(model) { }
  UInt operator() (const Element & element) {
    return 0;
  }
private:
  std::string mesh_data;
  const SolidMechanicsModel & model;
};

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_SELECTOR_HH__ */
