/**
 * @file   material_damage.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Chambart <marion.chambart@epfl.ch>
 * @date   Thu Jul 29 15:00:59 2010
 *
 * @brief  Material isotropic elastic
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
 *  if(key == "Yd") { sstr >> Yd; }
  else if(key == "Sd") { sstr >> Sd; }
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_MAZARS_HH__
#define __AKANTU_MATERIAL_MAZARS_HH__

__BEGIN_AKANTU__

/**
 * Material Mazars
 *
 * parameters in the material files :
 *   - rho  : density (default: 0)
 *   - E    : Young's modulus (default: 0)
 *   - nu   : Poisson's ratio (default: 1/2)
 *   - K0   : Damage threshold
 *   - At   : Parameter damage traction 1
 *   - Bt   : Parameter damage traction 2
 *   - Ac   : Parameter damage compression 1
 *   - Bc   : Parameter damage compression 2
 *   - beta : Parameter for shear
 */
class MaterialMazars : public virtual Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialMazars(Model & model, const ID & id = "");

  virtual ~MaterialMazars() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  virtual bool setParam(const std::string & key, const std::string & value,
			const ID & id);

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);



  /// Compute the tangent stiffness matrix for implicit for a given type
  void computeTangentStiffness(__attribute__ ((unused)) const ElementType & type,
			       __attribute__ ((unused)) Vector<double> & tangent_matrix,
			       __attribute__ ((unused)) GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };

  /// compute the stable time step for an element of size h
  Real getStableTimeStep(Real h, const Element & element = ElementNull) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;
  
protected:
  /// constitutive law for a given quadrature point
  __aka_inline__ void computeStress(Real * F, Real * sigma,Real & damage, Real & Ehat);
  
  __aka_inline__ void computeDamageAndStress( Real *F, Real * sigma,Real & damage, Real & Ehat);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Damage, damage, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the young modulus
  Real E;
  /// Poisson coefficient
  Real nu;

  /// damage threshold
  Real K0;
  ///parameter damage traction 1
  Real At ;
  ///parameter damage traction 2
  Real Bt ;
  ///parameter damage compression 1
  Real Ac ;
  ///parameter damage compression 2
  Real Bc ;
  ///parameter for shear
  Real beta ;

  /// damage internal variable
  ByElementTypeReal damage;
};

/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "material_mazars_inline_impl.cc"
#endif

/* -------------------------------------------------------------------------- */
/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MaterialMazars & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_MAZARS_HH__ */
