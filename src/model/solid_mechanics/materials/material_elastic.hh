/**
 * @file   material_elastic.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 29 15:00:59 2010
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

#ifndef __AKANTU_MATERIAL_ELASTIC_HH__
#define __AKANTU_MATERIAL_ELASTIC_HH__

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
class MaterialElastic : public virtual Material {
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

  virtual bool setParam(const std::string & key, const std::string & value,
		const ID & id);

  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  void computeTangentStiffness(const ElementType & el_type,
			       Vector<Real> & tangent_matrix,
			       GhostType ghost_type = _not_ghost);

  /// compute the p-wave speed in the material
  virtual Real getPushWaveSpeed();

  /// compute the s-wave speed in the material
  virtual Real getShearWaveSpeed();

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(types::Matrix & grad_u,
				  types::Matrix & sigma);

  // /// compute the tangent stiffness matrix for an element
  void computeTangentStiffnessOnQuad(types::Matrix & tangent);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the stable time step
  inline Real getStableTimeStep(Real h, const Element & element);

  AKANTU_GET_MACRO(E, E, Real);
  AKANTU_GET_MACRO(Nu, nu, Real);
  AKANTU_GET_MACRO(Mu, mu, Real);
  AKANTU_GET_MACRO(Lambda, lambda, Real);
  AKANTU_GET_MACRO(Kpa, kpa, Real);

  AKANTU_GET_MACRO(PlaneStress, plane_stress, bool);

  AKANTU_SET_MACRO(E, E, Real);
  AKANTU_SET_MACRO(Nu, nu, Real);

  virtual Real getParam(const ID & param) const;
  virtual void setParam(const ID & param, Real value);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// the young modulus
  Real E;

  /// Poisson coefficient
  Real nu;

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
