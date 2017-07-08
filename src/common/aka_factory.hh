/**
 * @file   aka_factory.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Thu Jul 06 2017
 *
 * @brief This is a generic factory
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
/* -------------------------------------------------------------------------- */
#include <map>
#include <memory>
#include <functional>
#include <string>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_FACTORY_HH__
#define __AKANTU_AKA_FACTORY_HH__

namespace akantu {

template <class Base, class... Args> class Factory {
  using allocator_t = std::function<std::unique_ptr<Base>(Args...)>;
private:
  Factory() = default;
public:
  Factory(const Factory &) = delete;
  Factory & operator=(const Factory &) = delete;

  static Factory & getInstance() {
    static Factory instance;
    return instance;
  }
  /* ------------------------------------------------------------------------ */
  bool registerAllocator(ID id, const allocator_t & allocator) {
    if (allocators.find(id) != allocators.end())
      AKANTU_EXCEPTION("The id " << id << " is already registered in the "
                       << debug::demangle(typeid(Base).name()) << " factory");
    allocators[id] = allocator;
    return true;
  }

  std::unique_ptr<Base> allocate(ID id, Args... args) {
    if (allocators.find(id) == allocators.end())
      AKANTU_EXCEPTION("The id  " << id << " is not registered in the "
                                  << debug::demangle(typeid(Base).name())
                                  << " factory.");
    return std::forward<std::unique_ptr<Base>>(allocators[id](args...));
  }

private:
  std::map<std::string, allocator_t> allocators;
};

} // namespace akantu

#endif /* __AKANTU_AKA_FACTORY_HH__ */
