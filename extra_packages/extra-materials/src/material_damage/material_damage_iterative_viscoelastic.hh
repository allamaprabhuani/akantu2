/**
 * @file  material_damage_iterative_viscoelastic.hh
 * @author  Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @date   Tue Nov 20 2018
 *
 * @brief  Header of the material iterative stiffness reduction viscoelastic
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
#include "material_iterative_stiffness_reduction.hh"
#include "material_viscoelastic_maxwell.hh"
/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_DAMAGE_ITERATIVE_VISCOELASTIC_HH__
#define __AKANTU_MATERIAL_DAMAGE_ITERATIVE_VISCOELASTIC_HH__

namespace akantu {

/**
 * Material damage iterative viscoelastic
 */

/* -------------------------------------------------------------------------- */
template <UInt dim>
class MaterialDamageIterativeViscoelastic
    : public MaterialIterativeStiffnessReductionIsotropic<
          dim, MaterialViscoelasticMaxwell> {
  using iterative_parent =
      MaterialIterativeStiffnessReductionIsotropic<dim,
                                                   MaterialViscoelasticMaxwell>;
  using viscous_grandparent = MaterialViscoelasticMaxwell<dim>;
  /* --------------------------------------------------------------------- */
  /* Constructors/Destructors                                              */
  /* --------------------------------------------------------------------- */
public:
  MaterialDamageIterativeViscoelastic(SolidMechanicsModel & model,
                                      const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// init the material
  void initMaterial() override;

  /// retract dissipated energy due to viscoelasticity
  /// and the mechanical work, save other internals
  void beforeSolveStep() override;

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  /// update the internal variables and compute energies dissipated due to
  /// damage and viscoelasticity
  void afterSolveStep(bool converged = true) override;

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(ElementType el_type, Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;

  /// update the last converged stress and strain arrays as well as current
  /// values of internal variables
  void updateIntVariables();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  using voigt_h = VoigtHelper<dim>;

  /// Values of viscous stress and stain at last converged step
  InternalField<Real> sigma_v_conv;
  InternalField<Real> epsilon_v_conv;

  /// values of converged grad U
  InternalField<Real> gradu_last;

  /// Dissipated energy
  InternalField<Real> dissipated_energy_damage;
};

} // namespace akantu

#endif /* __AKANTU_MATERIAL_DAMAGE_ITERATIVE_VISCOELASTIC_HH__ */
