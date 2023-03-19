/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
#include "material_orthotropic_damage.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_ITERATIVE_HH_
#define AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_ITERATIVE_HH_

namespace akantu {

/**
 * Material damage iterative
 *
 * parameters in the material files :
 *   - Sc
 */
template <Int spatial_dimension>
class MaterialOrthotropicDamageIterative
    : public MaterialOrthotropicDamage<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialOrthotropicDamageIterative(SolidMechanicsModel & model,
                                     const ID & id = "");

  virtual ~MaterialOrthotropicDamageIterative(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  ///  virtual void updateInternalParameters();

  virtual void computeAllStresses(GhostType ghost_type = _not_ghost);

  /// update internal field damage
  UInt updateDamage();

  /// update energies after damage has been updated
  virtual void updateEnergiesAfterDamage(ElementType el_type);

protected:
  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type,
                             GhostType ghost_type = _not_ghost);

  /// compute the equivalent stress on each Gauss point (i.e. the max prinicpal
  /// stress) and normalize it by the tensile strength
  void computeNormalizedEquivalentStress(const Array<Real> & grad_u,
                                         ElementType el_type,
                                         GhostType ghost_type = _not_ghost);

  /// find max normalized equivalent stress
  void findMaxNormalizedEquivalentStress(ElementType el_type,
                                         GhostType ghost_type = _not_ghost);

  inline void computeDamageAndStressOnQuad(Matrix<Real> & sigma,
                                           Matrix<Real> & one_minus_D,
                                           Matrix<Real> & root_one_minus_D,
                                           Matrix<Real> & damage,
                                           Matrix<Real> & first_term,
                                           Matrix<Real> & third_term);
  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get max normalized equivalent stress
  AKANTU_GET_MACRO(NormMaxEquivalentStress, norm_max_equivalent_stress, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// resistance to damage
  RandomInternalField<Real> Sc;

  /// internal field to store equivalent stress on each Gauss point
  InternalField<Real> equivalent_stress;

  /// internal field to store the direction of the current damage frame
  InternalField<Real> stress_dir;

  /// damage increment
  Real prescribed_dam;

  /// maximum equivalent stress
  Real norm_max_equivalent_stress;

  /// define damage threshold at which damage will be set to 1
  Real dam_threshold;

  /// quadrature point with highest equivalent Stress
  IntegrationPoint q_max;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_orthotropic_damage_iterative_inline_impl.hh"

} // namespace akantu

#endif /* AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_ITERATIVE_HH_ */
