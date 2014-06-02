/**
 * @file   aka_safe_enum.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Nov 28 17:15:04 2012
 *
 * @brief  Safe enums type (see More C++ Idioms/Type Safe Enum on Wikibooks
 * http://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Type_Safe_Enum)
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

#ifndef __AKANTU_AKA_SAFE_ENUM_HH__
#define __AKANTU_AKA_SAFE_ENUM_HH__

__BEGIN_AKANTU__
/// Safe enumerated type
template<typename def, typename inner = typename def::type>
class safe_enum : public def {
  typedef typename def::type type;
public:
  safe_enum(type v) : val(v) {}
  inner underlying() const { return val; }

  bool operator == (const safe_enum & s) const { return this->val == s.val; }
  bool operator != (const safe_enum & s) const { return this->val != s.val; }
  bool operator <  (const safe_enum & s) const { return this->val <  s.val; }
  bool operator <= (const safe_enum & s) const { return this->val <= s.val; }
  bool operator >  (const safe_enum & s) const { return this->val >  s.val; }
  bool operator >= (const safe_enum & s) const { return this->val >= s.val; }

  operator inner() { return val; }; 

public:
  // Works only if enumerations are contiguous.
  class iterator {
  public:
    iterator(type v) : it(v) { }
    void operator++() { ++it; }
    safe_enum operator*() { return static_cast<type>(it); }
    bool operator!=(iterator const & it) { return it.it != this->it; }
  private:
    int it;
  };

  static iterator begin() {
    return def::_begin_;
  }

  static iterator end() {
    return def::_end_;
  }

protected:
  inner val;
};

__END_AKANTU__

#endif /* __AKANTU_AKA_SAFE_ENUM_HH__ */
