/**
 * @file   material_cohesive_linear_friction.hh
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Tue Jul 29 2014
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * random stress definition for extrinsic type
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef __AKANTU_MATERIAL_COHESIVE_LINEAR_FRICTION_HH__
#define __AKANTU_MATERIAL_COHESIVE_LINEAR_FRICTION_HH__

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/**
 * Cohesive material linear with friction force
 *
 * parameters in the material files :
 *   - mu   : friction coefficient
 *   - penalty_for_friction : Penalty parameter for the friction behavior
 */
template<UInt spatial_dimension>
class MaterialCohesiveLinearFriction : public MaterialCohesiveLinear<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  typedef MaterialCohesiveLinear<spatial_dimension> MaterialParent;
public:

  MaterialCohesiveLinearFriction(SolidMechanicsModel & model, const ID & id = "");

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

  /// check delta_max for cohesive elements in case of no convergence
  /// in the solveStep (only for extrinsic-implicit)
  virtual void checkDeltaMax(GhostType ghost_type = _not_ghost);


  /// compute tangent stiffness matrix
  virtual void computeTangentTraction(const ElementType & el_type,
				      Array<Real> & tangent_matrix,
				      const Array<Real> & normal,
				      GhostType ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// maximum value of the friction coefficient
  Real mu_max;

  /// penalty parameter for the friction law
  Real friction_penalty;

  /// history parameter for the friction law
  CohesiveInternalField<Real> residual_sliding;

  /// friction force
  CohesiveInternalField<Real> friction_force;
};

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_COHESIVE_LINEAR_FRICTION_HH__ */
