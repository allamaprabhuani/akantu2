/**
 * @file   material_cohesive.hh
 *
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Wed Feb 22 16:31:20 2012
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
  typedef FEMTemplate<IntegratorGauss,
		      ShapeLagrange, _ek_cohesive>         MyFEMCohesiveType;
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
  virtual void resizeCohesiveArrays();

  /// compute tractions (including normals and openings)
  void computeTraction(GhostType ghost_type = _not_ghost);

  /// assemble residual
  void assembleResidual(GhostType ghost_type = _not_ghost);

  /// compute reversible and total energies by element
  void computeEnergies();

  /// check stress for cohesive elements' insertion
  virtual void checkInsertion(const Array<Real> & facet_stress,
			      const Mesh & mesh_facets,
			      ByElementTypeArray<bool> & facet_insertion);

  /// interpolate   stress  on   given   positions  for   each  element   (empty
  /// implemantation to avoid the generic call to be done on cohesive elements)
  virtual void interpolateStress(__attribute__((unused)) const ElementType type,
				 __attribute__((unused)) Array<Real> & result) { };

  virtual void computeAllStresses(__attribute__((unused)) GhostType ghost_type = _not_ghost) { };

  /// generate random sigma_c distributions
  void generateRandomDistribution(Array<Real> & sigma_lim);

protected:
  virtual void computeTangentTraction(__attribute__((unused)) const ElementType & el_type,
				      __attribute__((unused)) Array<Real> & tangent_matrix,
				      __attribute__((unused)) const Array<Real> & normal,
				      __attribute__((unused)) GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }


  void computeNormal(const Array<Real> & position,
		     Array<Real> & normal,
		     ElementType type,
		     GhostType ghost_type);

  void computeOpening(const Array<Real> & displacement,
		      Array<Real> & normal,
		      ElementType type,
		      GhostType ghost_type);

  template<ElementType type>
  void computeNormal(const Array<Real> & position,
		     Array<Real> & normal,
		     GhostType ghost_type);

  /// assemble stiffness
  void assembleStiffnessMatrix(GhostType ghost_type);

  /// constitutive law
  virtual void computeTraction(const Array<Real> & normal,
			       ElementType el_type,
			       GhostType ghost_type = _not_ghost) = 0;

  /// compute stress norms on quadrature points for each facet for stress check
  virtual void computeStressNorms(__attribute__((unused)) const Array<Real> & facet_stress,
				  __attribute__((unused)) Array<Real> & stress_check,
				  __attribute__((unused)) ElementType type_facet) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };


  /// parallelism functions
  inline UInt getNbDataForElements(const Array<Element> & elements,
				   SynchronizationTag tag) const;

  inline void packElementData(CommunicationBuffer & buffer,
			      const Array<Element> & elements,
			      SynchronizationTag tag) const;

  inline void unpackElementData(CommunicationBuffer & buffer,
				const Array<Element> & elements,
				SynchronizationTag tag);

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


  /// return the energy (identified by id) for the provided element
  virtual Real getEnergy(std::string energy_id, ElementType type, UInt index) {
    return Material::getEnergy(energy_id, type, index);
  }


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

  /// random generator
  RandomGenerator<Real> * random_generator;

  /// vector to store stresses on facets for element insertions
  Array<Real> sigma_insertion;

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
