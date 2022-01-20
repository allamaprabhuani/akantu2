/**
 * @file   material_drucker_prager.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Apr 06 2021
 *
 * @brief  Specialization of the material class for isotropic
 * plasticity with Drucker-Pruger yield criterion
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_voigthelper.hh"
#include "material_plastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_DRUCKER_PRAGER_HH__
#define __AKANTU_MATERIAL_DRUCKER_PRAGER_HH__

namespace akantu {

/**
 * Material plastic with a Drucker-pruger yield criterion
 */

template <Int dim> class MaterialDruckerPrager : public MaterialPlastic<dim> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  using Parent = MaterialPlastic<dim>;

public:
  MaterialDruckerPrager(SolidMechanicsModel & model, const ID & id = "");
  MaterialDruckerPrager(SolidMechanicsModel & model, Int spatial_dimension,
                        const Mesh & mesh, FEEngine & fe_engine,
                        const ID & id = "");

protected:
  using voigt_h = VoigtHelper<dim>;

  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

protected:
  template <class Args> inline void computeStressOnQuad(Args && args);
  inline void computeTangentModuliOnQuad(
      Matrix<Real> & tangent, const Matrix<Real> & grad_u,
      const Matrix<Real> & previous_grad_u, const Matrix<Real> & sigma_tensor,
      const Matrix<Real> & previous_sigma_tensor) const;

  inline Real computeYieldFunction(const Matrix<Real> & sigma);
  inline Real computeYieldStress(const Matrix<Real> & sigma);

  /// rcompute the alpha and k parameters
  void updateInternalParameters() override;

public:
  // closet point projection method to compute stress state on the
  // yield surface
  inline void computeGradientAndPlasticMultplier(
      const Matrix<Real, dim, dim> & sigma_tr, Real & plastic_multiplier_guess,
      Vector<Real, dim> & gradient_f,
      Vector<Real, dim> & delta_inelastic_strain, Int max_iterations = 100,
      Real tolerance = 1e-10);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  // Internal friction angle of the material
  Real phi;

  // Compressive strength of the material
  Real fc;

  // modified friction angle for Drucker-Prager
  Real alpha;

  // modified compressive strength for Drucker-Prager
  Real k;

  // radial return mapping
  bool radial_return_mapping;
};

} // namespace akantu

#include "material_drucker_prager_inline_impl.hh"

#endif /*__AKANTU_MATERIAL_DRUCKER_PRAGER_HH__  */
