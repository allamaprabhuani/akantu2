/**
 * @file   material_brittle.hh
 *
 * @author Josué Aranda <josue.arandaruiz@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Wed Feb 12 11:09:36 2014
 *
 * @brief  Brittle damage law
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
#include "material_damage.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_BRITTLE_HH__
#define __AKANTU_MATERIAL_BRITTLE_HH__

__BEGIN_AKANTU__

/**
 * Material brittle
 *
 * parameters in the material files :
 *   - S_0      : Critical stress at low strain rate (default: 157e6)
 *   - E_0      : Low strain rate threshold (default: 27e3)
 *   - A,B,C,D  : Fitting parameters for the critical stress at high strain rates 
 *                (default: 1.622e-11, -1.3274e-6, 3.6544e-2, -181.38)
 */
template<UInt spatial_dimension>
class MaterialBrittle : public MaterialDamage<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialBrittle(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialBrittle() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  virtual void updateInternalParameters();

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

protected:
  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(Matrix<Real> & grad_u,
                                  Matrix<Real> & grad_v,
				  Matrix<Real> & sigma,
				  Real & dam,
				  Real & sigma_equivalent);

  inline void computeDamageAndStressOnQuad(Matrix<Real> & sigma,
					   Real & dam,
					   Real & sigma_c,
					   Real & fracture_stress);

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:

  inline virtual UInt getNbDataForElements(const Array<Element> & elements,
					   SynchronizationTag tag) const;

  inline virtual void packElementData(CommunicationBuffer & buffer,
				      const Array<Element> & elements,
				      SynchronizationTag tag) const;

  inline virtual void unpackElementData(CommunicationBuffer & buffer,
					const Array<Element> & elements,
					SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// strain rate arrays ordered by element types
  InternalField<Real> strain_rate_brittle;

  //polynome constants for critical stress value
  Real A;
  Real B;
  Real C;
  Real D;

  //minimum strain rate
  Real E_0;

  //Critical stress at low strain rates
  Real S_0;

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_brittle_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_brittle_HH__ */
