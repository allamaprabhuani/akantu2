/**
 * @file   material_cohesive.hh
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Tue Feb  7 17:50:23 2012
 *
 * @brief  Specialization of the material class for cohesive elements
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
#include "fem_template.hh"
#include "aka_common.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_COHESIVE_HH__
#define __AKANTU_MATERIAL_COHESIVE_HH__

/* -------------------------------------------------------------------------- */
namespace akantu {
  class SolidMechanicsModelCohesive;
}

__BEGIN_AKANTU__

class MaterialCohesive : public Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef FEMTemplate< IntegratorCohesive<IntegratorGauss>,
		       ShapeCohesive<ShapeLagrange> >         MyFEMCohesiveType;
public:

  MaterialCohesive(SolidMechanicsModel& model, const ID & id = "");
  virtual ~MaterialCohesive();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize the material computed parameter
  virtual void initMaterial();

  /// resize vectors for new cohesive elements
  virtual void resizeCohesiveVectors();

  /// compute the residual for this material
  virtual void updateResidual(GhostType ghost_type = _not_ghost);

  /// check stress for cohesive elements' insertion
  virtual void checkInsertion(const Vector<Real> & facet_stress,
			      Vector<UInt> & facet_insertion);

  /// interpolate   stress  on   given   positions  for   each  element   (empty
  /// implemantation to avoid the generic call to be done on cohesive elements)
  virtual void interpolateStress(__attribute__((unused)) const ElementType type,
				 __attribute__((unused)) Vector<Real> & result) { };

protected:
  virtual void computeTangentTraction(__attribute__((unused)) const ElementType & el_type,
				      __attribute__((unused)) Vector<Real> & tangent_matrix,
				      __attribute__((unused)) const Vector<Real> & normal,
				      __attribute__((unused)) GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }


  void computeNormal(const Vector<Real> & position,
		     Vector<Real> & normal,
		     ElementType type,
		     GhostType ghost_type);

  void computeOpening(const Vector<Real> & displacement,
		      Vector<Real> & normal,
		      ElementType type,
		      GhostType ghost_type);

  template<ElementType type>
  void computeNormal(const Vector<Real> & position,
		     Vector<Real> & normal,
		     GhostType ghost_type);


  /// assemble residual
  void assembleResidual(GhostType ghost_type = _not_ghost);

  /// assemble stiffness
  void assembleStiffnessMatrix(GhostType ghost_type);

  /// compute tractions (including normals and openings)
  void computeTraction(GhostType ghost_type = _not_ghost);

  /// constitutive law
  virtual void computeTraction(const Vector<Real> & normal,
			       ElementType el_type,
			       GhostType ghost_type = _not_ghost) = 0;

  /// compute reversible and total energies by element
  void computeEnergies();

  /// compute stress norms on quadrature points for each facet for stress check
  virtual void computeStressNorms(__attribute__((unused)) const Vector<Real> & facet_stress,
				  __attribute__((unused)) Vector<Real> & stress_check) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// get the opening
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Opening, opening, Real);

  /// get the traction
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Traction, tractions, Real);

  /// compute reversible energy
  Real getReversibleEnergy();

  /// compute dissipated energy
  Real getDissipatedEnergy();

  /// get energy
  virtual Real getEnergy(std::string type);

  /// get damage
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Damage, damage, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// reversible energy by quadrature point
  ByElementTypeReal reversible_energy;

  /// total energy by quadrature point
  ByElementTypeReal total_energy;

  /// traction in all elements and quadrature points (previous time step)
  ByElementTypeReal tractions_old;

  /// opening in all elements and quadrature points (previous time step)
  ByElementTypeReal opening_old;

protected:

  /// traction in all elements and quadrature points
  ByElementTypeReal tractions;

  /// opening in all elements and quadrature points
  ByElementTypeReal opening;

  /// Link to the cohesive fem object in the model
  MyFEMCohesiveType * fem_cohesive;

  /// critical stress
  Real sigma_c;

  /// randomness factor
  Real rand;

  /// vector to store stresses on facets for element insertions
  Vector<Real> sigma_insertion;

  /// maximum displacement
  ByElementTypeReal delta_max;

  /// damage
  ByElementTypeReal damage;

  /// pointer to the solid mechanics model for cohesive elements
  SolidMechanicsModelCohesive * model;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_cohesive_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_COHESIVE_HH__ */
