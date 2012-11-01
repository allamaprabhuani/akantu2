/**
 * @file   material_marigo.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Feb 02 11:09:36 2012
 *
 * @brief  Marigo damage law
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

#ifndef __AKANTU_MATERIAL_MARIGO_HH__
#define __AKANTU_MATERIAL_MARIGO_HH__

__BEGIN_AKANTU__

/**
 * Material marigo
 *
 * parameters in the material files :
 *   - Yd  : (default: 50)
 *   - Sd  : (default: 5000)
 *   - Ydrandomness  : (default:0)
 */
template<UInt spatial_dimension>
class MaterialMarigo : public MaterialDamage<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialMarigo(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialMarigo() {};

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
  inline void computeStressOnQuad(types::RMatrix & grad_u,
				  types::RMatrix & sigma,
				  Real & dam,
				  Real & Y,
				  Real & Ydq);

  inline void computeDamageAndStressOnQuad(types::RMatrix & sigma,
					   Real & dam,
					   Real & Y,
					   Real & Ydq);

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:

  inline virtual UInt getNbDataForElements(const Vector<Element> & elements,
					   SynchronizationTag tag) const;

  inline virtual void packElementData(CommunicationBuffer & buffer,
				      const Vector<Element> & elements,
				      SynchronizationTag tag) const;

  inline virtual void unpackElementData(CommunicationBuffer & buffer,
					const Vector<Element> & elements,
					SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// resistance to damage
  Real Yd;

  /// damage threshold
  Real Sd;

  /// randomness on Yd
  Real Yd_randomness;

  /// critical epsilon when the material is considered as broken
  Real epsilon_c;

  Real Yc;

  /// Yd random internal variable
  ByElementTypeReal Yd_rand;

  bool damage_in_y;
  bool yc_limit;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_marigo_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_MARIGO_HH__ */
