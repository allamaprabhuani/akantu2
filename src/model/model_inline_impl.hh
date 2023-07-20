/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "model.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MODEL_INLINE_IMPL_HH_
#define AKANTU_MODEL_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename FEEngineClass>
inline FEEngineClass & Model::getFEEngineClassBoundary(const ID & name) {
  auto tmp_name = name.empty() ? default_fem : name;

  auto has_fem_boundary = hasFEEngineBoundary(tmp_name);
  auto boundary_name = tmp_name + ":boundary";

  if (not has_fem_boundary) {
    AKANTU_DEBUG_INFO("Creating FEEngine boundary " << tmp_name);

    auto & fe_engine = getFEEngine(tmp_name);

    auto spatial_dimension = fe_engine.getElementDimension();
    registerFEEngineObject<FEEngineClass>(boundary_name, fe_engine.getMesh(),
                                          spatial_dimension - 1);

    auto & fem_boundary = getFEEngineClass<FEEngineClass>(boundary_name);
    fem_boundary.computeNormalsOnIntegrationPoints(_not_ghost);
    fem_boundary.computeNormalsOnIntegrationPoints(_ghost);
    return fem_boundary;
  }

  return getFEEngineClass<FEEngineClass>(boundary_name);
}

/* -------------------------------------------------------------------------- */
template <typename FEEngineClass>
inline FEEngineClass & Model::getFEEngineClass(const ID & name) const {
  auto tmp_name = name.empty() ? default_fem : name;

  auto it = fems.find(tmp_name);
  if (it == fems.end()) {
    AKANTU_EXCEPTION("The FEEngine " << tmp_name << " is not registered");
  }

  return aka::as_type<FEEngineClass>(*(it->second));
}

/* -------------------------------------------------------------------------- */
inline void Model::unRegisterFEEngineObject(const ID & name) {
  auto it = fems.find(name);
  if (it == fems.end()) {
    AKANTU_EXCEPTION("FEEngine object with name " << name << " was not found");
  }

  fems.erase(it);
  if (not fems.empty() and default_fem == name) {
    default_fem = fems.begin()->first;
  }
}

/* -------------------------------------------------------------------------- */
template <typename FEEngineClass>
inline void Model::registerFEEngineObject(const ID & name, Mesh & mesh,
                                          Int spatial_dimension,
                                          bool do_not_precompute) {
  if (fems.empty()) {
    default_fem = name;
  }

  auto it = fems.find(name);
  if (it != fems.end()) {
    AKANTU_EXCEPTION("FEEngine object with name " << name
                                                  << " was already created");
  }

  fems[name] = std::make_unique<FEEngineClass>(
      mesh, spatial_dimension, id + ":fem:" + name, do_not_precompute);
}

/* -------------------------------------------------------------------------- */
inline FEEngine & Model::getFEEngine(const ID & name) const {
  ID tmp_name = (name.empty()) ? default_fem : name;

  auto it = fems.find(tmp_name);

  if (it == fems.end()) {
    AKANTU_EXCEPTION("The FEEngine " << tmp_name << " is not registered");
  }
  return *(it->second);
}

/* -------------------------------------------------------------------------- */
inline FEEngine & Model::getFEEngineBoundary(const ID & name) {
  ID tmp_name = (name.empty()) ? default_fem : name;

  auto it = fems.find(tmp_name + ":boundary");
  if (it == fems.end()) {
    AKANTU_EXCEPTION("The FEEngine boundary  " << tmp_name
                                               << " is not registered");
  }
  AKANTU_DEBUG_ASSERT(it->second != nullptr, "The FEEngine boundary "
                                                 << tmp_name
                                                 << " was not created");
  return *(it->second);
}

/* -------------------------------------------------------------------------- */
inline bool Model::hasFEEngineBoundary(const ID & name) {
  ID tmp_name = (name.empty()) ? default_fem : name;
  auto it = fems.find(tmp_name + ":boundary");
  return (it != fems.end());
}

/* -------------------------------------------------------------------------- */
template <typename T>
void Model::allocNodalField(std::unique_ptr<Array<T>> & array, Int nb_component,
                            const ID & name) const {
  if (array) {
    return;
  }

  auto nb_nodes = mesh.getNbNodes();
  array =
      std::make_unique<Array<T>>(nb_nodes, nb_component, T(), id + ":" + name);
}

/* -------------------------------------------------------------------------- */
inline Int Model::getNbIntegrationPoints(const Array<Element> & elements,
                                         const ID & fem_id) const {
  Int nb_quad = 0;
  for (auto && el : elements) {
    nb_quad +=
        getFEEngine(fem_id).getNbIntegrationPoints(el.type, el.ghost_type);
  }
  return nb_quad;
}

/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_MODEL_INLINE_IMPL_HH_ */
