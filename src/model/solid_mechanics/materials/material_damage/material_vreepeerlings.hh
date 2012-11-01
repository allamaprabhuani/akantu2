/**
 * @file   material_vreepeerlings.hh
 *
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date   Fri Feb 24 14:27:15 2012
 *
 * @brief  Specialization of the material class for the VreePeerlings material
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
#include "aka_common.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_VREEPEERLINGS_HH__
#define __AKANTU_MATERIAL_VREEPEERLINGS_HH__

__BEGIN_AKANTU__

/**
 * Material vreepeerlings
 *
 * parameters in the material files :
 *   - Kapa0i  : (default: 0.0001) Initial threshold (of the equivalent strain) for the initial step
 *   - Kapa0  : (default: 0.0001) Initial threshold (of the equivalent strain)
 *   - Alpha  : (default: 0.99) Fitting parameter (must be close to 1 to do tend to 0 the stress in the damaged element)
 *   - Beta   : (default: 300) This parameter determines the rate at which the damage grows 
 *   - Kct    : (default: 1) Ratio between compressive and tensile strength
 *   - Kapa0_randomness  : (default:0) Kapa random internal variable
 */
template<UInt spatial_dimension>
class MaterialVreePeerlings : public MaterialDamage<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialVreePeerlings(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialVreePeerlings() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

protected:
  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(types::RMatrix & F,
				  types::RMatrix & sigma,
				  Real & dam,
				  Real & Equistrain_rate,
				  Real & Equistrain,
				  Real & Kapaq,
				  Real dt,
				  types::RMatrix & strain_rate_vrpgls,
				  Real & crit_strain);

  inline void computeDamageAndStressOnQuad(types::RMatrix & sigma,
					   Real & dam,
					   Real & Equistrain_rate,
					   Real & Equistrain,
					   Real & Kapaq,
					   Real dt,
					   Real & crit_strain);

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Initial threshold (of the equivalent strain) (used in the initial step)
  Real Kapa0i;

  /// Initial threshold (of the equivalent strain)
  Real Kapa0;

  /// Fitting parameter (must be close to 1 to do tend to 0 the stress in the damaged element)
  Real Alpha;

  /// This parameter determines the rate at which the damage grows 
  Real Beta;

  /// Ratio between compressive and tensile strength
  Real Kct;

  /// randomness on Kapa0
  Real Kapa0_randomness;

  /// Kapa random internal variable
  ByElementTypeReal Kapa;

  /// strain_rate_vreepeerlings
  ByElementTypeReal strain_rate_vreepeerlings;

 /// strain_critical_vreepeerlings
  ByElementTypeReal critical_strain;

  /// Booleen to check the first step
  bool firststep;

  /// counter for initial step
  Int countforinitialstep;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_vreepeerlings_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_VREEPEERLINGS_HH__ */
