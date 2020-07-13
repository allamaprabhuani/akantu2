/**
 * @file   material_damage_iterative.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @date   Thu Feb 18 15:25:05 2016
 *
 * @brief  Damage material with constant stiffness reduction
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
#include "material_damage_iterative.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_DAMAGE_ITERATIVE_ORTHOTROPIC_HH__
#define __AKANTU_MATERIAL_DAMAGE_ITERATIVE_ORTHOTROPIC_HH__

namespace akantu {

/**
 * Material damage iterative orthotropic
 *
 */

template <UInt spatial_dimension>
class MaterialDamageIterativeOrthotropic
    : public MaterialDamageIterative<spatial_dimension,
                                     MaterialElasticOrthotropicHeterogeneous> {
  using parent =
      MaterialDamageIterative<spatial_dimension,
                              MaterialElasticOrthotropicHeterogeneous>;
  using OrthotropicParent =
      MaterialElasticOrthotropicHeterogeneous<spatial_dimension>;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialDamageIterativeOrthotropic(SolidMechanicsModel & model,
                                     const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial() override;

  void computeStress(const ElementType el_type, GhostType ghost_type) override;

  void computeTangentModuli(const ElementType & el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type) override;

protected:
  // recover damaged stiffness only in the direction where compression
  inline void reduceInternalParameters(Matrix<Real> & sigma, Real & dam,
                                       Real & _E1, Real & _E2, Real & _E3,
                                       Real & _nu12, Real & _nu13, Real & _nu23,
                                       Real & _G12, Real & _G13, Real & _G23,
                                       Matrix<Real> & _dir_vecs,
                                       Real & nb_flicks);
  // recover damaged stiffness based on the asymptotic function with
  // deformation normal to the crack
  inline void updateElasticModuli(Matrix<Real> & sigma, Matrix<Real> & grad_u,
                                  Real & dam, Real & _E1, Real & _E2,
                                  Real & _E3, Real & _nu12, Real & _nu13,
                                  Real & _nu23, Real & _G12, Real & _G13,
                                  Real & _G23, Matrix<Real> & _dir_vecs,
                                  Real & nb_flicks, bool & in_tension);

  // compute update moduli, compute C_prime and finally C
  inline void computeC(const ElementType el_type, GhostType ghost_type);

  // reset flickering counters before solve step
  void beforeSolveStep() override;
  /* ------------------------------------------------------------------------
   */
  /* DataAccessor inherited members */
  /* ------------------------------------------------------------------------
   */
  /* ------------------------------------------------------------------------
   */
  /* Accessors */
  /* ------------------------------------------------------------------------
   */
public:
  /* ------------------------------------------------------------------------
   */
  /* Class Members */
  /* ------------------------------------------------------------------------
   */
protected:
  using voigt_h = VoigtHelper<spatial_dimension>;

  /// number of state changes (kept as Real to be able to dump it)
  InternalField<Real> nb_state_changes;

  /// number of state changes
  InternalField<Real> damage_prev_iteration;

  /// stores previous stress state of an element
  InternalField<bool> in_tension;

  /// max allowed nb of state changes
  Real max_state_changes_allowed;

  /// flag to trigger contact
  bool contact;

  /// flag to apply stiffness in all directions
  bool iso_damage;

  /// flag to indicate material that it's in the loading test
  bool loading_test{false};

};

} // namespace akantu

#include "material_damage_iterative_orthotropic_inline_impl.cc"

#endif /* __AKANTU_MATERIAL_DAMAGE_ITERATIVE__HH__ */
