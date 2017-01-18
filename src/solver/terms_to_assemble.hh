/**
 * @file   terms_to_assemble.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Dec 20 2016
 *
 * @brief List of terms to assemble to a matrix
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
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TERMS_TO_ASSEMBLE_HH__
#define __AKANTU_TERMS_TO_ASSEMBLE_HH__

namespace akantu {

class TermsToAssemble {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  TermsToAssemble() {}
  virtual ~TermsToAssemble() {}

  class TermToAssemble {
  public:
    TermToAssemble(UInt i, UInt j) : _i(i), _j(j) {}
    inline void operator+=(Real val) { this->val = val; }
    inline operator Real() const  { return val; }
    inline UInt i() const { return _i; }
    inline UInt j() const { return _j; }
  private:
    UInt _i, _j;
    Real val;
  };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline TermToAssemble & operator()(UInt i, UInt j) {
    terms.emplace_back(i, j);
    return terms.back();
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
private:
  typedef std::vector<TermToAssemble> TermsContainer;
public:
  typedef TermsContainer::const_iterator const_terms_iterator;

  const_terms_iterator begin() const { return terms.begin(); }
  const_terms_iterator end() const { return terms.end(); }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  TermsContainer terms;
};

} // akantu

#endif /* __AKANTU_TERMS_TO_ASSEMBLE_HH__ */
