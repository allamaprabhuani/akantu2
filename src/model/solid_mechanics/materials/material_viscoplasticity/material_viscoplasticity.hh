/**
 * @file   material_plasticityinc.hh
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 *
 * @date   Tue Jul 09 18:15:37 20130
 *
 * @brief  Specialization of the material class for isotropic finite deformation linear hardening plasticityviscoplastic (small deformation)
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
#include "material_elastic.hh"
#include "aka_voigthelper.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_VISCOPLASTICITY_HH__
#define __AKANTU_MATERIAL_VISCOPLASTICITY_HH__

__BEGIN_AKANTU__

/**
 * Material plastic isotropic
 *
 * parameters in the material files :
 *   - h : Hardening parameter (default: 0)
 *   - sigmay : Yield stress
 *   - rate : Rate sensitivity
 *   - edot0 : Reference strain rate
 *
 *   - ts: Time step
 */


template <UInt spatial_dimension>
class MaterialViscoPlasticity : public MaterialElastic<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialViscoPlasticity(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialViscoPlasticity() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(const ElementType & el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost);

protected:
  inline void computeStressOnQuad(Matrix<Real> & grad_u,
				  Matrix<Real> & grad_delta_u,
				  Matrix<Real> & sigma,
				  Matrix<Real> & inelas_strain,
				  Real & iso_hardening);

  void computeTangentModuliOnQuad(Matrix<Real> & tangent,
				  Matrix<Real> & grad_delta_u,
				  Matrix<Real> & sigma_tensor,
				  Matrix<Real> & previous_sigma_tensor,
				  Real & iso_hardening);
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// Yield stresss
  Real sigmay;

  /// Hardening modulus
  Real h;

  /// Rate sensitivity component (rate)
  Real rate;

  /// Reference strain rate (edot0)
  Real edot0;

  /// Plane stress or plane strain
  bool plane_stress;

  /// Isotropic hardening (r)
  InternalField<Real> iso_hardening;

  /// Time step (ts)
  Real ts;

};
/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "material_viscoplasticity_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_VISCOPLASTICITY_HH__ */
