/**
 * @file   material_damage.hh
 *
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Mar 15 16:06:20 2011
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

#ifndef __AKANTU_MATERIAL_DAMAGE_HH__
#define __AKANTU_MATERIAL_DAMAGE_HH__

__BEGIN_AKANTU__
template<UInt spatial_dimension, template<UInt> class Parent = MaterialElastic>
class MaterialDamage : public Parent<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialDamage(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialDamage() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  virtual void computeAllStresses(GhostType ghost_type);

protected:
  /// update the dissipated energy, must be called after the stress have been computed
  void updateDissipatedEnergy(GhostType ghost_type);


  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// give the dissipated energy for the time step
  Real getDissipatedEnergy() const;

  virtual Real getEnergy(std::string type);
  virtual Real getEnergy(std::string energy_id, ElementType type, UInt index) {
    return Parent<spatial_dimension>::getEnergy(energy_id, type, index);
  };

  AKANTU_GET_MACRO_NOT_CONST(Damage, damage, ByElementTypeReal &);
  AKANTU_GET_MACRO(Damage, damage, const ByElementTypeReal &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// damage internal variable
  ByElementTypeReal damage;

  /// dissipated energy
  ByElementTypeReal dissipated_energy;

  /// previous strain used to compute the dissipated energy
  ByElementTypeReal strain_prev;

  /// previous strain used to compute the dissipated energy
  ByElementTypeReal stress_prev;

  /// contain the current value of @f$ \int_0^{\epsilon}\sigma(\omega)d\omega @f$ the dissipated energy
  ByElementTypeReal int_sigma;
};


__END_AKANTU__

#endif /* __AKANTU_MATERIAL_DAMAGE_HH__ */
