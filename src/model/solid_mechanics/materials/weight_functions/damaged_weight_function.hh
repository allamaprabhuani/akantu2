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
#include "base_weight_function.hh"
/* -------------------------------------------------------------------------- */
#ifndef AKANTU_DAMAGED_WEIGHT_FUNCTION_HH_
#define AKANTU_DAMAGED_WEIGHT_FUNCTION_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
/*  Damage weight function                                                    */
/* -------------------------------------------------------------------------- */
class DamagedWeightFunction : public BaseWeightFunction {
public:
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  DamagedWeightFunction(NonLocalManager & manager)
      : BaseWeightFunction(manager, "damaged") {
    this->init();
  }

  /* ------------------------------------------------------------------------ */
  /* Base Weight Function inherited methods                                   */
  /* ------------------------------------------------------------------------ */
  /// set the pointers of internals to the right flattend version
  void init() override;

  inline Real operator()(Real r,
                         const __attribute__((unused)) IntegrationPoint & q1,
                         const IntegrationPoint & q2);

private:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

  /// internal pointer to the current damage vector
  ElementTypeMapReal * damage{nullptr};
};

} // namespace akantu

#include "damaged_weight_function_inline_impl.hh"

#endif /* AKANTU_DAMAGED_WEIGHT_FUNCTION_HH_ */
