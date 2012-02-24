/**
 * @file   material_vreeperlings.hh
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 * @date   Fri Feb 17 14:00:00 2012
 *
 * @brief  Specialization of the material class for the VreePerlings material
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

#ifndef __AKANTU_MATERIAL_VREEPERLINGS_HH__
#define __AKANTU_MATERIAL_VREEPERLINGS_HH__

__BEGIN_AKANTU__

/**
 * Material vreeperlings
 *
 * parameters in the material files :
 *   - Kapa0  : (default: 50) Initial threshold (of the equivalent strain)
 *   - Alpha  : (default: 0.99) Fitting parameter (must be close to 1 to do tend to 0 the stress in the damaged element)
 *   - Beta   : (default: 300) This parameter determines the rate at which the damage grows 
 *   - Kct    : (default: 1) Ratio between compressive and tensile strength
 *   - Kapa0_randomness  : (default:0) Kapa random internal variable
 */
class MaterialVreePerlings : public MaterialDamage {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialVreePerlings(Model & model, const ID & id = "");

  virtual ~MaterialVreePerlings() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  bool setParam(const std::string & key, const std::string & value,
		const ID & id);

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// constitutive law for a given quadrature point
  __aka_inline__ void computeStress(Real * F, Real * sigma, Real & vreeperlings, Real & Equistrain, Real & Kapaq);

  __aka_inline__ void computeDamageAndStress(Real * sigma, Real & dam, Real & Equistrain, Real & Kapaq );

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:

  __aka_inline__ virtual UInt getNbDataToPack(const Element & element,
 				      SynchronizationTag tag);

  __aka_inline__ virtual UInt getNbDataToUnpack(const Element & element,
 					SynchronizationTag tag);

  __aka_inline__ virtual void packData(CommunicationBuffer & buffer,
 			       const Element & element,
 			       SynchronizationTag tag);

  __aka_inline__ virtual void unpackData(CommunicationBuffer & buffer,
                                 const Element & element,
                                 SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

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

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "material_vreeperlings_inline_impl.cc"
#endif

/* -------------------------------------------------------------------------- */
/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MaterialVreePerlings & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_VREEPERLINGS_HH__ */
