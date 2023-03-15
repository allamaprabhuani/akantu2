/**
 * Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Regularisation>
NTNFricLawLinearSlipWeakeningNoHealing<Regularisation>::
    NTNFricLawLinearSlipWeakeningNoHealing(NTNBaseContact & contact,
                                           const ID & id)
    : NTNFricLawLinearSlipWeakening<Regularisation>(contact, id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearSlipWeakeningNoHealing<
    Regularisation>::computeFrictionCoefficient() {
  AKANTU_DEBUG_IN();

  // get arrays
  const SynchronizedArray<Real> & slip = this->internalGetCumulativeSlip();

  auto nb_contact_nodes = this->contact.getNbContactNodes();
  for (Int n = 0; n < nb_contact_nodes; ++n) {
    if (slip(n) >= this->d_c(n)) {
      this->mu(n) = this->mu_k(n);
    } else {
      // mu = mu_k + (1 - slip / Dc) * (mu_s - mu_k)
      this->mu(n) = this->mu_k(n) + (1 - (slip(n) / this->d_c(n))) *
                                        (this->mu_s(n) - this->mu_k(n));
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawLinearSlipWeakeningNoHealing<Regularisation>::printself(
    std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT) {
    ;
  }

  stream << space << "NTNFricLawLinearSlipWeakeningNoHealing [" << std::endl;
  NTNFricLawLinearSlipWeakening<Regularisation>::printself(stream, ++indent);
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
