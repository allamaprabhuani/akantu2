/**
 * @file   material_cohesive_linear.hh
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Mon Feb 20 12:00:34 2012
 *
 * @brief Linear irreversible cohesive law of mixed mode loading with
 * random stress definition for extrinsic type
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
#include "material_cohesive.hh"
#include "aka_common.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_COHESIVE_LINEAR_EXTRINSIC_HH__
#define __AKANTU_MATERIAL_COHESIVE_LINEAR_EXTRINSIC_HH__

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/**
 * Cohesive material linear damage for extrinsic case
 *
 * parameters in the material files :
 *   - sigma_c   : critical stress sigma_c  (default: 0)
 *   - beta      : weighting parameter for sliding and normal opening (default: 0)
 *   - G_cI      : fracture energy for mode I (default: 0)
 *   - G_cII     : fracture energy for mode II (default: 0)
 *   - rand      : randomness factor (default: 0)
 *   - penalty   : stiffness in compression to prevent penetration
 */
template<UInt spatial_dimension>
class MaterialCohesiveLinearExtrinsic : public MaterialCohesive {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialCohesiveLinearExtrinsic(SolidMechanicsModel & model, const ID & id = "");
  virtual ~MaterialCohesiveLinearExtrinsic();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize the material computed parameter
  virtual void initMaterial();

  /// resize vectors for new cohesive elements
  virtual void resizeCohesiveVectors();

  /// compute stress norms on quadrature points for each facet for stress check
  virtual void computeStressNorms(const Vector<Real> & facet_stress,
				  Vector<Real> & stress_check);

protected:

  /// constitutive law
  void computeTraction(const Vector<Real> & normal,
		       ElementType el_type,
		       GhostType ghost_type = _not_ghost);

  /// compute effective stress norm for insertion check
  inline Real computeEffectiveNorm(const types::Matrix & stress,
				   const types::RVector & normal,
				   const types::RVector & tangent);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// critical effective stress
  ByElementTypeReal sigma_c_eff;

  /// beta parameter
  Real beta;

  /// mode I fracture energy
  Real G_cI;

  /// mode II fracture energy
  Real G_cII;

  /// kappa parameter
  Real kappa;

  /// penalty coefficient
  Real penalty;

  /// critical displacement
  ByElementTypeReal delta_c;

  /// vector to temporarily store the normal stress for the norm
  types::RVector normal_stress;

  /// vector to temporarily store the tangential stress for the norm
  types::RVector tangential_stress;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "material_cohesive_linear_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_COHESIVE_LINEAR_EXTRINSIC_HH__ */
