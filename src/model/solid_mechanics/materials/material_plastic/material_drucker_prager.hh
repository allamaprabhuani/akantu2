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
 * Material plastic with a Drucker-Prager yield criterion
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
  void computeTangentModuli(ElementType el_type,
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

  inline Real computeYieldStress(const Matrix<Real> & sigma);

  inline void computeDeviatoricStress(const Matrix<Real> & sigma,
				      Matrix<Real> & sigma_dev);

  inline bool aboveThresholdStress(const Matrix<Real> &,
				   Real & , Real & );

  inline void projectionOnThreshold(Matrix<Real> &, Vector<Real> &,
				    Real &, Real &);

  /// compute the alpha and k parameters
  void updateInternalParameters() override;


 
  /* ------------------------------------------------------------------------ */
  /* Methods for projection to the yield surface                              */
  /* ------------------------------------------------------------------------ */  
  /// compute the objective function to minimize
  virtual inline Real computeObjectiveFunction(const Matrix<Real> & sigma_guess,
					       const Matrix<Real> & sigma_trial,
					       Real & plastic_multiplier_guess,
					       Vector<Real> & gradient_f,
					       Vector<Real> & delta_inelastic_strain,
					       Real & yield_function);

  /// compute the jacobian of the objective function
  virtual inline void computeJacobian(const Matrix<Real> &, Vector<Real> &);

  /// compute the hessian of the objective function
  virtual inline void computeHessian(const Matrix<Real> &, Matrix<Real> &);


  
public:
  // closet point projection method to compute stress state on the
  // yield surface
  inline void computeGradientAndPlasticMultplier(
      const Matrix<Real> & sigma_tr, Real & plastic_multiplier_guess,      
      Vector<Real> & gradient_f, Vector<Real> & delta_inelastic_strain,
      UInt max_iterations = 100, Real tolerance = 1e-8);

public:
    /// get modified compressive strength
  AKANTU_GET_MACRO(K, k, Real);

  /// get modified friction angle
  AKANTU_GET_MACRO(Alpha, alpha, Real);

  
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

  // maximum number of iterations for projection
  UInt max_iterations;

  // tolerance error for projection
  Real tolerance;
  
  // maximum deviatoric for capped Drucker-Prager
  Real deviatoric_max;

  // minimum hydrostatic for capped Drucker-Prager
  Real hydrostatic_min;

  // capped
  bool is_capped;
};
  
} // namespace akantu


#include "material_drucker_prager_inline_impl.hh"


#endif /*__AKANTU_MATERIAL_DRUCKER_PRAGER_HH__  */
