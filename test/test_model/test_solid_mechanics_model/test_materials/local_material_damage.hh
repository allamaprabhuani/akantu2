/**
 * @file   local_material_damage.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Nov 26 00:17:56 2010
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

#ifndef __AKANTU_LOCAL_MATERIAL_DAMAGE_HH__
#define __AKANTU_LOCAL_MATERIAL_DAMAGE_HH__

__BEGIN_AKANTU__

class LocalMaterialDamage : public Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  LocalMaterialDamage(SolidMechanicsModel & model, const ID & id = "");

  virtual ~LocalMaterialDamage() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(types::RMatrix & grad_u,
				  types::RMatrix & sigma,
				  Real & damage);

  /// compute tangent stiffness
  virtual void computeTangentStiffness(__attribute__ ((unused)) const ElementType & el_type,
				       __attribute__ ((unused)) Vector<Real> & tangent_matrix,
				       __attribute__ ((unused)) GhostType ghost_type = _not_ghost) {};

  /// compute the potential energy for all elements
  void computePotentialEnergy(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the potential energy for on element
  inline void computePotentialEnergyOnQuad(types::RMatrix & grad_u,
					   types::RMatrix & sigma,
					   Real & epot);

  /// compute the celerity of wave in the material
  inline Real celerity();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the stable time step
  inline Real getStableTimeStep(Real h, const Element & element);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Damage, damage, Real);
private:

  /// the young modulus
  Real E;

  /// Poisson coefficient
  Real nu;

  /// First Lamé coefficient
  Real lambda;

  /// Second Lamé coefficient (shear modulus)
  Real mu;

  /// resistance to damage
  Real Yd;

  /// damage threshold
  Real Sd;

  /// Bulk modulus
  Real kpa;

  /// damage internal variable
  ByElementTypeReal damage;

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "local_material_damage_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_LOCAL_MATERIAL_DAMAGE_HH__ */
