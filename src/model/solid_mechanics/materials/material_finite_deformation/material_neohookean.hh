/**
 * @file   material_neohookean.hh
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Wed Aug 04 10:58:42 2010
 *
 * @brief  Material isotropic elastic
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
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_NEOHOOKEAN_HH__
#define __AKANTU_MATERIAL_NEOHOOKEAN_HH__

__BEGIN_AKANTU__

/**
 * Material elastic isotropic
 *
 * parameters in the material files :
 *   - rho : density (default: 0)
 *   - E   : Young's modulus (default: 0)
 *   - nu  : Poisson's ratio (default: 1/2)
 *   - Plane_Stress : if 0: plane strain, else: plane stress (default: 0)
 */
template<UInt spatial_dimension>
class MaterialNeohookean : public virtual Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialNeohookean(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialNeohookean() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void initMaterial();

  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the elastic potential energy
  virtual void computePotentialEnergy(ElementType el_type,
				      GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(const ElementType & el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost);

  /// compute the p-wave speed in the material
  virtual Real getPushWaveSpeed(const Element & element) const;

  /// compute the s-wave speed in the material
  virtual Real getShearWaveSpeed(const Element & element) const;

protected:

  /// constitutive law for a given quadrature point
  inline void computePiolaKirchhoffOnQuad(const Matrix<Real> & E,
                                          Matrix<Real> & S);

  /// constitutive law for a given quadrature point
  inline void computeDeltaStressOnQuad(const Matrix<Real> & grad_u, const Matrix<Real> & grad_delta_u,
        Matrix<Real> & delta_S);

  inline void computeStressOnQuad(Matrix<Real> & grad_u,
                                  Matrix<Real> & sigma);

  /// constitutive law for a given quadrature point
  //inline void updateStressOnQuad(const Matrix<Real> & sigma,
  //                              Matrix<Real> & cauchy_sigma);

  /// compute the potential energy for a quadrature point
  inline void computePotentialEnergyOnQuad(const Matrix<Real> & grad_u,
                                           Real & epot);

  /// compute the tangent stiffness matrix for an element
  void computeTangentModuliOnQuad(Matrix<Real> & tangent, Matrix<Real> & grad_u);

  /// recompute the lame coefficient if E or nu changes
  virtual void updateInternalParameters();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// the young modulus
  Real E;

  /// Poisson coefficient
  Real nu;

  /// Neohookean 2
  /// the young modulus at infinite strain
  Real E1;
  /// the young modulus  decay parameter
  Real alpha;
  /// Second Lamé coefficient (shear modulus) at infinite strain
  Real mu1;
  /// First Lamé coefficient at infinite strain
  Real lambda1;

  /// First Lamé coefficient
  Real lambda;

  /// Second Lamé coefficient (shear modulus)
  Real mu;

  /// Bulk modulus
  Real kpa;

  /// Plane stress or plane strain
  bool plane_stress;

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_neohookean_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_NEOHOOKEAN_HH__ */
