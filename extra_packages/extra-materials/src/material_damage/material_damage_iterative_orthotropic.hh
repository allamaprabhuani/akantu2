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
  inline void computeOrthotropicStress(
      Matrix<Real> & sigma, Matrix<Real> & grad_u, Real & dam, Real & _E1,
      Real & _E2, Real & _E3, Real & _nu12, Real & _nu13, Real & _nu23,
      Real & _G12, Real & _G13, Real & _G23, Matrix<Real> & _Cprime,
      Matrix<Real> & _C, Vector<Real> & _eigC, Matrix<Real> & _dir_vecs,
      UInt & nb_flicks);

  // reset flickering counters before solve step
  void beforeSolveStep();
  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  using voigt_h = VoigtHelper<spatial_dimension>;

  /// initial elastic modulus = E1 = E2 = E3
  Real E;

  /// number of times element changed from tension to compression
  InternalField<UInt> nb_state_changes;

  /// flag responsible to fix stiffness if element is flickering
  bool fix_flickering_elements;

  /// max allowed nb of state changes
  UInt max_state_changes_allowed;
};

} // namespace akantu

#include "material_damage_iterative_orthotropic_inline_impl.cc"

#endif /* __AKANTU_MATERIAL_DAMAGE_ITERATIVE__HH__ */
