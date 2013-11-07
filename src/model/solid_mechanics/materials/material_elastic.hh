/**
 * @file   material_elastic.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
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
#include "material_thermal.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_ELASTIC_HH__
#define __AKANTU_MATERIAL_ELASTIC_HH__

__BEGIN_AKANTU__

/**
 * Material elastic isotropic
 *
 * parameters in the material files :
 *   - E   : Young's modulus (default: 0)
 *   - nu  : Poisson's ratio (default: 1/2)
 *   - Plane_Stress : if 0: plane strain, else: plane stress (default: 0)
 */
template<UInt spatial_dimension>
class MaterialElastic : public MaterialThermal<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialElastic(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialElastic() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void initMaterial();

  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(const ElementType & el_type,
			    Array<Real> & tangent_matrix,
			    GhostType ghost_type = _not_ghost);

  /// compute the p-wave speed in the material
  virtual Real getPushWaveSpeed() const;

  /// compute the s-wave speed in the material
  virtual Real getShearWaveSpeed() const;

protected:
  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(const Matrix<Real> & grad_u,
				  Matrix<Real> & sigma, Real sigma_th = 0);

  /// compute the tangent stiffness matrix for an element
  void computeTangentModuliOnQuad(Matrix<Real> & tangent);

  /// recompute the lame coefficient if E or nu changes
  virtual void updateInternalParameters();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the stable time step
  inline Real getStableTimeStep(Real h, const Element & element);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

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

#include "material_elastic_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_ELASTIC_HH__ */
