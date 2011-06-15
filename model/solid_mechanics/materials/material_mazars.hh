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
 *
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

class MaterialMazars : public Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialMazars(Model & model, const MaterialID & id = "");

  virtual ~MaterialMazars() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  void setParam(const std::string & key, const std::string & value,
		const MaterialID & id);

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// constitutive law for a given quadrature point
  inline void computeStress(Real * F, Real * sigma,Real & damage);

  /// compute the potential energy for all elements
  void computePotentialEnergy(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the potential energy for on element
  inline void computePotentialEnergy(Real * F, Real * sigma, Real * epot);

  /// Compute the tangent stiffness matrix for implicit for a given type
  void computeTangentStiffness(__attribute__ ((unused)) const ElementType & type,
			       __attribute__ ((unused)) Vector<double> & tangent_matrix,
			       __attribute__ ((unused)) GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };
    
  /// compute the celerity of wave in the material
  inline Real celerity();

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the stable time step
  inline Real getStableTimeStep(Real h, const Element & element);

  /// return damage value
  ByElementTypeReal & getDamage(){return damage;};

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Damage, damage, const Vector<Real> &);
private:

  /// the young modulus
  Real E;

  /// Poisson coefficient
  Real nu;

  /// First Lamé coefficient
  Real lambda;

  /// Second Lamé coefficient (shear modulus)
  Real mu;


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


  /// Bulk modulus
  Real kpa;

  /// damage internal variable
  ByElementTypeReal damage;
  ByElementTypeReal ghost_damage;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_mazars_inline_impl.cc"

/* -------------------------------------------------------------------------- */
/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MaterialMazars & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_MAZARS_HH__ */
