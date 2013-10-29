/**
 * @file   material_neohookean.hh
 *
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Jun 15 10:45:12 2011
 *
 * @brief  Material Neo Hookean
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
 *   - Plain_Stress : if 0: plain strain, else: plain stress (default: 0)
 */
template<UInt spatial_dimension>
class MaterialNeohookean : public MaterialElastic<spatial_dimension> {
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

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(const ElementType & el_type,
			    Array<Real> & tangent_matrix,
			    GhostType ghost_type = _not_ghost);

  /// compute the potential energy for all elements
  virtual void computePotentialEnergy(ElementType el_type, GhostType ghost_type = _not_ghost);

private:
  /// compute the celerity of wave in the material
  Real celerity(const Element & element);

  /// compute the potential energy for on element
  inline void computePotentialEnergyOnQuad(Matrix<Real> & grad_u,
					   Real & epot);
protected:
  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(Matrix<Real> & grad_u,
				  Matrix<Real> & sigma);

  // /// compute the tangent stiffness matrix for an element
  void computeTangentModuliOnQuad(Matrix<Real> & grad_u,
				  Matrix<Real> & tangent);

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

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_neohookean_inline_impl.cc"


__END_AKANTU__

#endif /* __AKANTU_MATERIAL_NEOHOOKEAN_HH__ */
