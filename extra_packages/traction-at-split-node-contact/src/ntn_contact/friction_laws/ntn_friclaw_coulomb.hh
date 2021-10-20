/**
 * @file   ntn_friclaw_coulomb.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Sep 29 2020
 *
 * @brief  coulomb friction with \mu_s = \mu_k (constant)
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef AST_NTN_FRICLAW_COULOMB_HH_
#define AST_NTN_FRICLAW_COULOMB_HH_

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_fricreg_no_regularisation.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Regularisation = NTNFricRegNoRegularisation>
class NTNFricLawCoulomb : public Regularisation {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NTNFricLawCoulomb(NTNBaseContact & contact, const ID & id = "coulomb");
  ~NTNFricLawCoulomb() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// register synchronizedarrays for sync
  void registerSynchronizedArray(SynchronizedArrayBase & array) override;

  /// dump restart file
  void dumpRestart(const std::string & file_name) const override;

  /// read restart file
  void readRestart(const std::string & file_name) override;

  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

protected:
  /// compute frictional strength according to friction law
  void computeFrictionalStrength() override;

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  void addDumpFieldToDumper(const std::string & dumper_name,
                            const std::string & field_id) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // friction coefficient
  SynchronizedArray<Real> mu;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
template <class Regularisation>
inline std::ostream &
operator<<(std::ostream & stream,
           const NTNFricLawCoulomb<Regularisation> & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "ntn_friclaw_coulomb_tmpl.hh"

#endif /* AST_NTN_FRICLAW_COULOMB_HH_ */
