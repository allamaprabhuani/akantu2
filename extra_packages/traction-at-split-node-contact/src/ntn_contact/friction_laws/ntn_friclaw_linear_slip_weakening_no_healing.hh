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
#ifndef AST_NTN_FRICLAW_LINEAR_SLIP_WEAKENING_NO_HEALING_HH_
#define AST_NTN_FRICLAW_LINEAR_SLIP_WEAKENING_NO_HEALING_HH_

/* -------------------------------------------------------------------------- */
#include "ntn_friclaw_linear_slip_weakening.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Regularisation = NTNFricRegNoRegularisation>
class NTNFricLawLinearSlipWeakeningNoHealing
    : public NTNFricLawLinearSlipWeakening<Regularisation> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NTNFricLawLinearSlipWeakeningNoHealing(
      NTNBaseContact & contact,
      const ID & id = "linear_slip_weakening_no_healing");
  ~NTNFricLawLinearSlipWeakeningNoHealing() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

protected:
  /// computes the friction coefficient as a function of slip
  void computeFrictionCoefficient() override;

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
template <class Regularisation>
inline std::ostream & operator<<(
    std::ostream & stream,
    const NTNFricLawLinearSlipWeakeningNoHealing<Regularisation> & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "ntn_friclaw_linear_slip_weakening_no_healing_tmpl.hh"

#endif /* AST_NTN_FRICLAW_LINEAR_SLIP_WEAKENING_NO_HEALING_HH_ */
