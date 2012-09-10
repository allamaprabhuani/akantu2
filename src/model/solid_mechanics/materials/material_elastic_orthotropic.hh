/**
 * @file   material_elastic_orthotropic.hh
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Thu Apr 12 10:57:52 2012
 *
 * @brief  Orthotropic elastic material
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_ELASTIC_ORTHOTROPIC_HH__
#define __AKANTU_MATERIAL_ELASTIC_ORTHOTROPIC_HH__

__BEGIN_AKANTU__

/**
 * Orthotropic elastic material
 * \todo extend this law for any material orientation
 *
 * parameters in the material files :
 *   - rho  : density (default: 0)
 *   - E1   : Young's modulus along x (default: 0)
 *   - E2   : Young's modulus along y (default: 0)
 *   - E3   : Young's modulus along z (default: 0)
 *   - nu12 : Poisson's ratio along xy (default: 0)
 *   - nu13 : Poisson's ratio along xz (default: 0)
 *   - nu23 : Poisson's ratio along yz (default: 0)
 *   - G12  : Shear modulus along xy (default: 0)
 *   - G13  : Shear modulus along xz (default: 0)
 *   - G23  : Shear modulus along yz (default: 0)
 *   - Plane_Stress : if 0: plane strain, else: plane stress (default: 0)
 */
template<UInt spatial_dimension>
class MaterialElasticOrthotropic : public virtual Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialElasticOrthotropic(SolidMechanicsModel & model, const ID & id = "");

  ~MaterialElasticOrthotropic();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void initMaterial();

  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(const ElementType & el_type,
			    Vector<Real> & tangent_matrix,
			    GhostType ghost_type = _not_ghost);

  /// compute the p-wave speed in the material
  virtual Real getPushWaveSpeed() const;

  /// compute the s-wave speed in the material
  virtual Real getShearWaveSpeed() const;

  virtual void updateInternalParameters();
protected:
  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(types::Matrix & grad_u,
				  types::Matrix & sigma);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the stable time step
  inline Real getStableTimeStep(Real h, const Element & element);

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

  /// the x young modulus
  Real E1;

  /// the y young modulus
  Real E2;

  /// the z young modulus
  Real E3;

  /// xy Poisson coefficient
  Real nu12;

  /// xz Poisson coefficient
  Real nu13;

  /// yz Poisson coefficient
  Real nu23;

  /// xy shear modulus
  Real G12;

  /// xz shear modulus
  Real G13;

  /// yz shear modulus
  Real G23;

  /// stiffness coefficients
  types::Matrix * S;

  /// Plane stress or plane strain
  bool plane_stress;

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_elastic_orthotropic_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_ELASTIC_ORTHOTROPIC_HH__ */
