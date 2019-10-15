/**
 * @file   material_elastic_orthotropic.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Enrico Milanese <enrico.milanese@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Feb 16 2018
 *
 * @brief  Orthotropic elastic material
 *
 * @section LICENSE
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material_elastic_linear_anisotropic_heterogeneous.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_ELASTIC_ORTHOTROPIC_HETEROGENEOUS_HH__
#define __AKANTU_MATERIAL_ELASTIC_ORTHOTROPIC_HETEROGENEOUS_HH__

namespace akantu {

/**
 * Orthotropic elastic material
 *
 * parameters in the material files :
 *   - n1   : direction of x-axis in material base, normalisation not necessary
 * (default: {1, 0, 0})
 *   - n2   : direction of y-axis in material base, normalisation not necessary
 * (default: {0, 1, 0})
 *   - n3   : direction of z-axis in material base, normalisation not necessary
 * (default: {0, 0, 1})
 *   - rho  : density (default: 0)
 *   - E1   : Young's modulus along n1 (default: 0)
 *   - E2   : Young's modulus along n2 (default: 0)
 *   - E3   : Young's modulus along n3 (default: 0)
 *   - nu12 : Poisson's ratio along 12 (default: 0)
 *   - nu13 : Poisson's ratio along 13 (default: 0)
 *   - nu23 : Poisson's ratio along 23 (default: 0)
 *   - G12  : Shear modulus along 12 (default: 0)
 *   - G13  : Shear modulus along 13 (default: 0)
 *   - G23  : Shear modulus along 23 (default: 0)
 */

template <UInt Dim>
class MaterialElasticOrthotropicHeterogeneous
    : public MaterialElasticLinearAnisotropicHeterogeneous<Dim> {
  using parent = MaterialElasticLinearAnisotropicHeterogeneous<Dim>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialElasticOrthotropicHeterogeneous(SolidMechanicsModel & model,
                                          const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial() override;

  void updateInternalParameters() override;
  void updateInternalParametersOnQuad(const Real & _E1, const Real & _E2,
                                      const Real & _E3, const Real & _nu12,
                                      const Real & _nu13, const Real & _nu23,
                                      const Real & _G12, const Real & _G13,
                                      const Real & _G23, Matrix<Real> & _Cprime,
                                      Matrix<Real> & _C, Vector<Real> & _eigC,
                                      Matrix<Real> & _dir_vecs);

  void
  computePotentialEnergyByElement(ElementType type, UInt index,
                                  Vector<Real> & epot_on_quad_points) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(E1, E1, Real);
  AKANTU_GET_MACRO(E2, E2, Real);
  AKANTU_GET_MACRO(E3, E3, Real);
  AKANTU_GET_MACRO(Nu12, nu12, Real);
  AKANTU_GET_MACRO(Nu13, nu13, Real);
  AKANTU_GET_MACRO(Nu23, nu23, Real);
  AKANTU_GET_MACRO(G12, G12, Real);
  AKANTU_GET_MACRO(G13, G13, Real);
  AKANTU_GET_MACRO(G23, G23, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  using voigt_h = VoigtHelper<Dim>;
  /// the n1 young modulus
  Real E1;

  /// the n2 young modulus
  Real E2;

  /// the n3 young modulus
  Real E3;

  /// 12 Poisson coefficient
  Real nu12;

  /// 13 Poisson coefficient
  Real nu13;

  /// 23 Poisson coefficient
  Real nu23;

  /// 12 shear modulus
  Real G12;

  /// 13 shear modulus
  Real G13;

  /// 23 shear modulus
  Real G23;

  InternalField<Real> E1_field;
  InternalField<Real> E2_field;
  InternalField<Real> E3_field;
  InternalField<Real> nu12_field;
  InternalField<Real> nu13_field;
  InternalField<Real> nu23_field;
  InternalField<Real> G12_field;
  InternalField<Real> G13_field;
  InternalField<Real> G23_field;

  bool are_internals_initialized;
};

} // namespace akantu

#endif /* __AKANTU_MATERIAL_ELASTIC_ORTHOTROPIC_HETEROGENEOUS_HH__ */
