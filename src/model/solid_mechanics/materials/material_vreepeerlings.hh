/**
 * @file   material_vreepeerlings.hh
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 * @date   Fri Feb 17 14:00:00 2012
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
 *   - Kapa0  : (default: 50) Initial threshold (of the equivalent strain)
 *   - Alpha  : (default: 0.99) Fitting parameter (must be close to 1 to do tend to 0 the stress in the damaged element)
 *   - Beta   : (default: 300) This parameter determines the rate at which the damage grows 
 *   - Kct    : (default: 1) Ratio between compressive and tensile strength
 *   - Kapa0_randomness  : (default:0) Kapa random internal variable
 */
class MaterialVreePeerlings : public MaterialDamage {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialVreePeerlings(Model & model, const ID & id = "");

  virtual ~MaterialVreePeerlings() {};

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
  inline void computeStress(Real * F, Real * sigma, Real & dam, Real & Equistrain, Real & Kapaq);

  inline void computeDamageAndStress(Real * sigma, Real & dam, Real & Equistrain, Real & Kapaq);

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:

  virtual UInt getNbDataToPack(const Element & element,
			       SynchronizationTag tag);

  virtual UInt getNbDataToUnpack(const Element & element,
				 SynchronizationTag tag);

  virtual void packData(CommunicationBuffer & buffer,
			const Element & element,
			SynchronizationTag tag);

  virtual void unpackData(CommunicationBuffer & buffer,
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

#include "material_vreepeerlings_inline_impl.cc"

/* -------------------------------------------------------------------------- */
/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MaterialVreePeerlings & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_VREEPEERLINGS_HH__ */
