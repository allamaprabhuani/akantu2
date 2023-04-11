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
#include "aka_common.hh"
#include "resolution_penalty.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_RESOLUTION_PENALTY_QUADRATIC_HH_
#define AKANTU_RESOLUTION_PENALTY_QUADRATIC_HH_

namespace akantu {

class ResolutionPenaltyQuadratic : public ResolutionPenalty {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
private:
  using Parent = ResolutionPenalty;

public:
  ResolutionPenaltyQuadratic(ContactMechanicsModel & model, const ID & id = "");

  ~ResolutionPenaltyQuadratic() override = default;

public:
  /// local computation of tangent moduli due to normal traction
  void computeNormalModuli(const ContactElement & element,
                           Matrix<Real> & stiffness) override;

protected:
  /// local computation of normal traction due to penetration
  Real computeNormalTraction(const Real & gap) const override;
};

} // namespace akantu

#endif /* AKANTU_RESOLUTION_PENALTY_QUADRATIC_HH_ */
