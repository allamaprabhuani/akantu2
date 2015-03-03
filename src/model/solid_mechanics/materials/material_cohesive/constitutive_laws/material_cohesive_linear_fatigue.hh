/**
 * @file   material_cohesive_linear_fatigue.hh
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Thu Feb 19 14:20:59 2015
 *
 * @brief Linear irreversible cohesive law with dissipative
 * unloading-reloading cycles
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

#include "material_cohesive_linear.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_COHESIVE_LINEAR_FATIGUE_HH__
#define __AKANTU_MATERIAL_COHESIVE_LINEAR_FATIGUE_HH__

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/**
 * Linear irreversible cohesive law with dissipative
 * unloading-reloading cycles
 *
 * This law uses two different stiffnesses during unloading and
 * reloading. The implementation is based on the article entitled "A
 * cohesive model for fatigue crack growth" by Nguyen, Repetto, Ortiz
 * and Radovitzky (2001). This law is identical to the
 * MaterialCohesiveLinear one except for the unloading-reloading
 * phase.
 *
 * input parameter:
 *
 * - delta_f : it must be greater than delta_c and it is inversely
 *      proportional to the dissipation in the unloading-reloading
 *      cycles (default: delta_c)
 */

template <UInt spatial_dimension>
class MaterialCohesiveLinearFatigue : public MaterialCohesiveLinear<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialCohesiveLinearFatigue(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize the material parameters
  virtual void initMaterial();

protected:

  /// constitutive law
  virtual void computeTraction(const Array<Real> & normal,
			       ElementType el_type,
			       GhostType ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// delta_f parameter
  Real delta_f;

  /// delta of the previous step
  CohesiveInternalField<Real> delta_prec;

  /// stiffness for reloading
  CohesiveInternalField<Real> K_plus;

  /// stiffness for unloading
  CohesiveInternalField<Real> K_minus;

  /// 1D traction in the cohesive law
  CohesiveInternalField<Real> T_1d;
};

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_COHESIVE_LINEAR_FATIGUE_HH__ */
