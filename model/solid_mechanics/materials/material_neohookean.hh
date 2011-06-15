/**
 * @file   material_neohookean.hh
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
class MaterialNeohookean : public Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialNeohookean(Model & model, const MaterialID & id = "");

  virtual ~MaterialNeohookean() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  void setParam(const std::string & key, const std::string & value,
		const MaterialID & id);

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  void computeTangentStiffness(const ElementType & el_type,
			       Vector<Real> & tangent_matrix,
			       GhostType ghost_type = _not_ghost);

  /// compute the potential energy for all elements
  void computePotentialEnergy(ElementType el_type, GhostType ghost_type = _not_ghost);


  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  /// compute the celerity of wave in the material
  Real celerity(const Element & element);

  /// compute the potential energy for on element
  inline void computePotentialEnergy(Real * F, Real * epot);
protected:
  /// constitutive law for a given quadrature point
  inline void computeStress(Real * F, Real * sigma);

  // /// compute the tangent stiffness matrix for an element type
  template<UInt dim>
  void computeTangentStiffnessByDim(akantu::ElementType, akantu::Vector<Real>& tangent_matrix, akantu::GhostType);

  // /// compute the tangent stiffness matrix for an element
  template<UInt dim>
  void computeTangentStiffness(Real * tangent);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the stable time step
  inline Real getStableTimeStep(Real h, const Element & element);
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

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

  /// Plain stress or plain strain
  bool plain_stress;

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_neohookean_inline_impl.cc"

/* -------------------------------------------------------------------------- */
/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MaterialNeohookean & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_NEOHOOKEAN_HH__ */
