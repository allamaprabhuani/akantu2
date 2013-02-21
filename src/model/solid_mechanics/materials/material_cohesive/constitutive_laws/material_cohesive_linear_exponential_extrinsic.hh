/**
 * @file   material_cohesive_linear_exponential_extrinsic.hh
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Thu May 24 10:46:59 2012
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
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
#ifndef __AKANTU_MATERIAL_COHESIVE_LINEAR_EXPONENTIAL_EXTRINSIC_HH__
#define __AKANTU_MATERIAL_COHESIVE_LINEAR_EXPONENTIAL_EXTRINSIC_HH__

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/**
 * Cohesive material linear damage for extrinsic case
 *
 * parameters in the material files :
 *   - sigma_c   : critical stress  (default: 0)
 *   - sigma_max : maximum stress  (default: 0)
 *   - delta_0   : delta 0  (default: 0)
 *   - beta      : weighting parameter for sliding and normal opening (default: 0)
 *   - G_cI      : fracture energy for mode I (default: 0)
 *   - G_cII     : fracture energy for mode II (default: 0)
 *   - rand      : randomness factor (default: 0)
 *   - gamma     : gamma factor (default: 0)
 *   - z_max     : z max (default: 0)
 */
template<UInt spatial_dimension>
class MaterialCohesiveLinearExponentialExtrinsic : public MaterialCohesive {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialCohesiveLinearExponentialExtrinsic(SolidMechanicsModel & model, const ID & id = "");
  virtual ~MaterialCohesiveLinearExponentialExtrinsic();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize the material computed parameter
  virtual void initMaterial();

  /// resize vectors for new cohesive elements
  virtual void resizeCohesiveVectors();

  /// compute effective stress norm for insertion check
  virtual Real computeEffectiveNorm(const types::RMatrix & stress,
				    const types::RVector & normal,
				    const types::RVector & tangent);

protected:

  /// constitutive law
  void computeTraction(const Vector<Real> & normal,
		       ElementType el_type,
		       GhostType ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// maximum stress
  Real sigma_max;

  /// delta 0
  Real delta_0;

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

  /// actual maximum sigma
  ByElementTypeReal sigma_actual;

  /// gamma factor
  Real gamma;

  /// z max
  Real z_max;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "material_cohesive_linear_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_COHESIVE_LINEAR_EXPONENTIAL_EXTRINSIC_HH__ */
