﻿/**
 * Copyright (©) 2012-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_cohesive_damage.hh"
//#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <numeric>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
MaterialCohesiveDamage::MaterialCohesiveDamage(SolidMechanicsModel & model,
                                                    const ID & id)
    : MaterialCohesive(model, id)
{
  AKANTU_DEBUG_IN();

  this->registerParam("k", k, Real(0.), _pat_parsable | _pat_readable,
                      "Cohesive stiffness parameter");

  this->registerParam("G_c", G_c, Real(0.), _pat_parsable | _pat_readable,
                      "Mode I fracture energy");
  a = 2.*k*G_c/(sigma_c*sigma_c);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesiveDamage::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialCohesive::initMaterial();

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
Real MaterialCohesiveDamage::stiffness(Real d)
{
    return k*(1./d-1.);
}
Real MaterialCohesiveDamage::augmented_stiffness(Real d)
{
    return k*1./d;
}
Real MaterialCohesiveDamage::augmented_compliance(Real d)
{
    return d/k;
}

}
