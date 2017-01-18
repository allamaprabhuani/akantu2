/**
 * @file   material_damage.hh
 *
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sat Jul 11 2015
 *
 * @brief  Material isotropic elastic
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_DAMAGE_HH__
#define __AKANTU_MATERIAL_DAMAGE_HH__

__BEGIN_AKANTU__
template <UInt spatial_dimension,
          template <UInt> class Parent = MaterialElastic>
class MaterialDamage : public Parent<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialDamage(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialDamage(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void initMaterial();

  /// compute the tangent stiffness matrix for an element type
  virtual void computeTangentModuli(const ElementType & el_type,
                                    Array<Real> & tangent_matrix,
                                    GhostType ghost_type = _not_ghost);
  virtual bool hasStiffnessMatrixChanged() { return true; }

protected:
  /// update the dissipated energy, must be called after the stress have been
  /// computed
  virtual void updateEnergies(ElementType el_type, GhostType ghost_type);

  /// compute the tangent stiffness matrix for a given quadrature point
  inline void computeTangentModuliOnQuad(Matrix<Real> & tangent, Real & dam);

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

  AKANTU_GET_MACRO_NOT_CONST(Damage, damage, ElementTypeMapArray<Real> &);
  AKANTU_GET_MACRO(Damage, damage, const ElementTypeMapArray<Real> &);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Damage, damage, Real)

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// damage internal variable
  InternalField<Real> damage;

  /// dissipated energy
  InternalField<Real> dissipated_energy;

  /// contain the current value of @f$ \int_0^{\epsilon}\sigma(\omega)d\omega
  /// @f$ the dissipated energy
  InternalField<Real> int_sigma;
};

__END_AKANTU__

#include "material_damage_tmpl.hh"

#endif /* __AKANTU_MATERIAL_DAMAGE_HH__ */
