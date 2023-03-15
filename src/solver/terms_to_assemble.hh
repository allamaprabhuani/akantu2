/**
 * @file   terms_to_assemble.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Oct 11 2017
 *
 * @brief  List of terms to assemble to a matrix
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
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_TERMS_TO_ASSEMBLE_HH_
#define AKANTU_TERMS_TO_ASSEMBLE_HH_

namespace akantu {

class TermsToAssemble {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  TermsToAssemble(const ID & dof_id_m, const ID & dof_id_n)
      : dof_id_m(dof_id_m), dof_id_n(dof_id_n) {}
  virtual ~TermsToAssemble() = default;

  class TermToAssemble {
  public:
    TermToAssemble(Idx i, Idx j) : _i(i), _j(j) {}
    inline TermToAssemble & operator=(Real val) {
      this->val = val;
      return *this;
    }
    inline TermToAssemble operator+=(Real val) {
      this->val += val;
      return *this;
    }
    inline operator Real() const { return val; }
    inline Idx i() const { return _i; }
    inline Idx j() const { return _j; }

  private:
    Idx _i, _j;
    Real val{0.};
  };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline TermToAssemble & operator()(Idx i, Idx j) {
    terms.emplace_back(i, j);
    return terms.back();
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
private:
  using TermsContainer = std::vector<TermToAssemble>;

public:
  using const_terms_iterator = TermsContainer::const_iterator;

  auto begin() const { return terms.begin(); }
  auto end() const { return terms.end(); }

  AKANTU_GET_MACRO(DOFIdM, dof_id_m, const ID &);
  AKANTU_GET_MACRO(DOFIdN, dof_id_n, const ID &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  TermsContainer terms;
  ID dof_id_m, dof_id_n;
};

} // namespace akantu

#endif /* AKANTU_TERMS_TO_ASSEMBLE_HH_ */
