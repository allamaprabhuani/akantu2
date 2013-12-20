/**
 * @file   material_plasticityinc.hh
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date   Tue Jul 09 18:15:37 20130
 *
 * @brief  Specialization of the material class for isotropic finite deformation linear hardening plasticity
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
#include "aka_voigthelper.hh"
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_PLASTICITYINC_HH__
#define __AKANTU_MATERIAL_PLASTICITYINC_HH__

__BEGIN_AKANTU__


/**
 * Material plastic isotropic
 *
 * parameters in the material files :
 *   - h : Hardening parameter (default: 0)
 *   - sigmay : Yield stress
 */
template <UInt spatial_dimension>
class MaterialPlasticityinc : public MaterialElastic<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialPlasticityinc(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialPlasticityinc() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// Compute the stresses and the plastic strain energy
  virtual void computeAllStresses(GhostType ghost_type = _not_ghost);

  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(const ElementType & el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost);

  virtual Real getPlasticEnergy();
  virtual Real getEnergy(std::string type);

  /// Compute the plastic energy
  virtual void updatePlasticEnergy(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// Compute the plastic energy increment
  virtual void computePlasticEnergyIncrement(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// Compute the true potential energy
  virtual void computePotentialEnergy(ElementType el_type, GhostType ghost_type);

protected:

  /// constitutive law for a given quadrature point
  //inline void computePiolaKirchhoffOnQuad(const Matrix<Real> & E,
  //Matrix<Real> & S);

  /// constitutive law for a given quadrature point
  //inline void computeDeltaStressOnQuad(const Matrix<Real> & grad_u, const Matrix<Real> & grad_delta_u,
  //      Matrix<Real> & delta_S);

  inline void computeStressOnQuad(Matrix<Real> & grad_u,
                                  Matrix<Real> & grad_delta_u,
                                  Matrix<Real> & sigma,
                                  Matrix<Real> & inelas_strain,
                                  Matrix<Real> & d_inelas_strain,
                                  Real & iso_hardening,
                                  Real sigma_th_cur,
                                  Real sigma_th_prev);

  inline void computeTangentModuliOnQuad(Matrix<Real> & tangent,
                                         Matrix<Real> & grad_delta_u,
                                         Matrix<Real> & sigma_tensor,
                                         Matrix<Real> & previous_sigma_tensor,
                                         Real & iso_hardening);

  inline void computePotentialEnergyOnQuad(Matrix<Real> & grad_u,
                                           Matrix<Real> & sigma,
                                           const Matrix<Real> & inelas_strain,
                                           Real & epot);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// Yield stresss
  Real sigmay;

  /// hardening modulus
  Real h;

  /// isotropic hardening, r
  InternalField<Real> iso_hardening;

  /// Plastic energy
  InternalField<Real> plastic_energy;

  /// Plastic energy increment
  InternalField<Real> d_plastic_energy;

  /// Plastic strain increment
  InternalField<Real> d_inelas_strain;
};
/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "material_plasticityinc_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_PLASTICITYINC_HH__ */
