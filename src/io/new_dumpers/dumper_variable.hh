/**
 * @file   dumper_variable.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Jun 04 2013
 * @date last modification: Sun Oct 19 2014
 *
 * @brief  template of variable
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_DUMPER_VARIABLE_HH__
#define __AKANTU_DUMPER_VARIABLE_HH__
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace dumper {
  /* ------------------------------------------------------------------------ */
  /// Variable interface
  class VariableBase {
  public:
    VariableBase() = default;
    virtual ~VariableBase() = default;
  };

  /* ------------------------------------------------------------------------ */
  template <typename T, bool is_scal = std::is_scalar<T>::value>
  class Variable : public VariableBase {
  public:
    explicit Variable(const T & t) : var(t) {}

    const T & operator[](UInt i) const { return var[i]; }

    UInt getDim() { return var.size(); }

  protected:
    const T & var;
  };

  /* ------------------------------------------------------------------------ */
  template <typename T> class Variable<T, true> : public VariableBase {
  public:
    explicit Variable(const T & t) : var(t) {}

    const T & operator[](UInt) const { return var; }
    UInt getDim() { return 1; }

  protected:
    const T & var;
  };

} // namespace dumper
} // namespace akantu

#endif /* __AKANTU_DUMPER_VARIABLE_HH__ */
