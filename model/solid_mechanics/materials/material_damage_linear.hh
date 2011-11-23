/**
 * @file   material_damage.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marion Chambart <marion.chambart@epfl.ch>
 * @date   Thu Jul 29 15:00:59 2010
 *
 * @brief  Material isotropic elastic + linear softening
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

#ifndef __AKANTU_MATERIAL_DAMAGE_LINEAR_HH__
#define __AKANTU_MATERIAL_DAMAGE_LINEAR_HH__

__BEGIN_AKANTU__

/**
 * Material liner damage
 *
 * parameters in the material files :
 *   - Sigc : (default: 1e5)
 *   - Gc  : (default: 2)
 */
class MaterialDamageLinear : public MaterialElastic {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialDamageLinear(Model & model, const ID & id = "");

  virtual ~MaterialDamageLinear() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  bool setParam(const std::string & key, const std::string & value,
		const ID & id);

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// Compute the tangent stiffness matrix for implicit for a given type
  void computeTangentStiffness(__attribute__ ((unused)) const ElementType & type,
			       __attribute__ ((unused)) Vector<double> & tangent_matrix,
			       __attribute__ ((unused)) GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };




  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// constitutive law for a given quadrature point
__aka_inline__ void computeStress(Real * F, Real * sigma, Real & damage, Real &K);

  virtual __aka_inline__ Real getStableTimeStep(Real h, const Element & element) {
    return MaterialElastic::getStableTimeStep(h, element);
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Damage, damage, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// kind of toughness 
  Real Gc;

  /// critical stress
  Real Sigc;

  /// damage internal variable
  ByElementTypeReal damage;
  ByElementTypeReal K; 

Real Epsmin, Epsmax;
};

/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "material_damage_linear_inline_impl.cc"
#endif

/* -------------------------------------------------------------------------- */
/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MaterialDamageLinear & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_DAMAGE_LINEAR_HH__ */
