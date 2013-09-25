/**
 * @file   material_elastic_orthotropic.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   Tue May 08 13:01:18 2012
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

#ifndef __AKANTU_MATERIAL_ELASTIC_LINEAR_ANISOTROPIC_HH__
#define __AKANTU_MATERIAL_ELASTIC_LINEAR_ANISOTROPIC_HH__

__BEGIN_AKANTU__

/**
 * General linear anisotropic elastic material
 * The only constraint on the elastic tensor is that it can be represented
 * as a symmetric 6x6 matrix (3D) or 3x3 (2D).
 *
 * parameters in the material files :
 *   - rho  : density (default: 0)
 *   - C_ij  : entry on the stiffness
 */
template<UInt Dim>
class MaterialElasticLinearAnisotropic : public virtual Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialElasticLinearAnisotropic(SolidMechanicsModel & model, const ID & id = "",
                                   bool symmetric = true);

  ~MaterialElasticLinearAnisotropic();

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



  virtual void updateInternalParameters();


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// compute stable timestep
  virtual Real getStableTimeStep(Real h, const Element & element);
  /// compute max wave celerity
  Real getCelerity() const;


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// stiffness coefficients
  Matrix<Real>  C;
  /// eigenvalues of stiffness tensor
  Vector<Real> eigC;
  const static VoigtHelper<Dim> voigt_h;
  bool symmetric;

};
__END_AKANTU__

#endif /* __AKANTU_MATERIAL_ELASTIC_LINEAR_ANISOTROPIC_HH__ */
