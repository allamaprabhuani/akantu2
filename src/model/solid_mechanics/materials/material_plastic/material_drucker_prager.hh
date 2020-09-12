/**
 * @file   material_drucker_pruger.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed Sep 09 2020
 * @date last modification: Wed Sep 09 2020
 *
 * @brief  Specialization of the material class for isotropic
 * plasticity with Drucker-Pruger yield criterion
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
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
  
template<UInt spatial_dimension>
class MaterialDruckerPrager
  : public MaterialPlastic<spatial_dimension>  {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialDruckerPrager(SolidMechanicsModel & model,
			const ID & id = "");
  MaterialDruckerPrager(SolidMechanicsModel & model, UInt dim,
			const Mesh & mesh, FEEngine & fe_engine,
			const ID & id = "");

protected:
  using voigt_h = VoigtHelper<spatial_dimension>;
  
  void initialize();

 
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(const ElementType & el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;

protected:
  /// Infinitesimal deformations
  inline void computeStressOnQuad(
      const Matrix<Real> & grad_u, const Matrix<Real> & previous_grad_u,
      Matrix<Real> & sigma, const Matrix<Real> & previous_sigma,
      Matrix<Real> & inelas_strain, const Matrix<Real> & previous_inelas_strain,
      const Real & sigma_th, const Real & previous_sigma_th);

  /// Finite deformations
  inline void computeStressOnQuad(
      const Matrix<Real> & grad_u, const Matrix<Real> & previous_grad_u,
      Matrix<Real> & sigma, const Matrix<Real> & previous_sigma,
      Matrix<Real> & inelas_strain, const Matrix<Real> & previous_inelas_strain,
      const Real & sigma_th, const Real & previous_sigma_th,
      const Matrix<Real> & F_tensor);

  inline void computeTangentModuliOnQuad(
      Matrix<Real> & tangent, const Matrix<Real> & grad_u,
      const Matrix<Real> & previous_grad_u, const Matrix<Real> & sigma_tensor,
      const Matrix<Real> & previous_sigma_tensor) const;

  inline Real computeYieldFunction(const Matrix<Real> & sigma);

  inline void computeDeviatoricStress(const Matrix<Real> & sigma,
				      Matrix<Real> & sigma_dev);

public:
  // closet point projection method to compute stress state on the
  // yield surface
  inline void computeGradientAndPlasticMultplier(
      const Matrix<Real> & sigma_tr, Real & plastic_multiplier_guess,      
      Vector<Real> & gradient_f, Vector<Real> & delta_inelastic_strain,
      UInt max_iterations = 100, Real tolerance = 1e-10);

  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  // friction angle
  Real phi;

  // compressive strength
  Real fc;

  //
  Real alpha;

  //
  Real k;
};
  
} // namespace akantu


#include "material_drucker_prager_inline_impl.hh"


#endif /*__AKANTU_MATERIAL_DRUCKER_PRAGER_HH__  */
