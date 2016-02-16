/**
 * @file   material_cohesive.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Jan 12 2016
 *
 * @brief  Specialization of the material class for cohesive elements
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "fe_engine_template.hh"
#include "aka_common.hh"
#include "cohesive_internal_field.hh"
#include "cohesive_element_inserter.hh"

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
  typedef FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_cohesive>
      MyFEEngineCohesiveType;

public:
  MaterialCohesive(SolidMechanicsModel & model, const ID & id = "");
  virtual ~MaterialCohesive();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material computed parameter
  virtual void initMaterial();

  /// compute tractions (including normals and openings)
  void computeTraction(GhostType ghost_type = _not_ghost);

  /// assemble residual
  void assembleInternalForces(GhostType ghost_type = _not_ghost);

  /// compute reversible and total energies by element
  void computeEnergies();

  /// check stress for cohesive elements' insertion, by default it
  /// also updates the cohesive elements' data
  virtual void checkInsertion(__attribute__((unused)) bool check_only = false) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// check delta_max for cohesive elements in case of no convergence
  /// in the solveStep (only for extrinsic-implicit)
  virtual void checkDeltaMax(__attribute__((unused))
                             GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// reset variables when convergence is reached (only for
  /// extrinsic-implicit)
  virtual void resetVariables(__attribute__((unused))
                              GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// interpolate   stress  on   given   positions  for   each  element   (empty
  /// implemantation to avoid the generic call to be done on cohesive elements)
  virtual void interpolateStress(__attribute__((unused)) const ElementType type,
                                 __attribute__((unused))
                                 Array<Real> & result){};

  /// compute the stresses
  virtual void computeAllStresses(__attribute__((unused))
                                  GhostType ghost_type = _not_ghost){};

  // add the facet to be handled by the material
  UInt addFacet(const Element & element);

protected:
  virtual void
  computeTangentTraction(__attribute__((unused)) const ElementType & el_type,
                         __attribute__((unused)) Array<Real> & tangent_matrix,
                         __attribute__((unused)) const Array<Real> & normal,
                         __attribute__((unused))
                         GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  /// compute the normal
  void computeNormal(const Array<Real> & position, Array<Real> & normal,
                     ElementType type, GhostType ghost_type);

  /// compute the opening
  void computeOpening(const Array<Real> & displacement, Array<Real> & normal,
                      ElementType type, GhostType ghost_type);

  template <ElementType type>
  void computeNormal(const Array<Real> & position, Array<Real> & normal,
                     GhostType ghost_type);

  /// assemble stiffness
  void assembleStiffnessMatrix(GhostType ghost_type);

  /// constitutive law
  virtual void computeTraction(const Array<Real> & normal, ElementType el_type,
                               GhostType ghost_type = _not_ghost) = 0;

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

  /// get damage
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Damage, damage, Real);

  /// get facet filter
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(FacetFilter, facet_filter, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(FacetFilter, facet_filter, UInt);
  AKANTU_GET_MACRO(FacetFilter, facet_filter,
                   const ElementTypeMapArray<UInt> &);
  // AKANTU_GET_MACRO(ElementFilter, element_filter, const
  // ElementTypeMapArray<UInt> &);

  /// compute reversible energy
  Real getReversibleEnergy();

  /// compute dissipated energy
  Real getDissipatedEnergy();

  /// compute contact energy
  Real getContactEnergy();

  /// get energy
  virtual Real getEnergy(std::string type);

  /// return the energy (identified by id) for the provided element
  virtual Real getEnergy(std::string energy_id, ElementType type, UInt index) {
    return Material::getEnergy(energy_id, type, index);
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// list of facets assigned to this material
  ElementTypeMapArray<UInt> facet_filter;

  /// Link to the cohesive fem object in the model
  MyFEEngineCohesiveType * fem_cohesive;

private:
  /// reversible energy by quadrature point
  CohesiveInternalField<Real> reversible_energy;

  /// total energy by quadrature point
  CohesiveInternalField<Real> total_energy;

protected:
  /// opening in all elements and quadrature points
  CohesiveInternalField<Real> opening;

  /// opening in all elements and quadrature points (previous time step)
  CohesiveInternalField<Real> opening_old;

  /// traction in all elements and quadrature points
  CohesiveInternalField<Real> tractions;

  /// traction in all elements and quadrature points (previous time step)
  CohesiveInternalField<Real> tractions_old;

  /// traction due to contact
  CohesiveInternalField<Real> contact_tractions;

  /// normal openings for contact tractions
  CohesiveInternalField<Real> contact_opening;

  /// maximum displacement
  CohesiveInternalField<Real> delta_max;

  /// tell if the previous delta_max state is needed (in iterative schemes)
  bool use_previous_delta_max;

  /// tell if the previous opening state is needed (in iterative schemes)
  bool use_previous_opening;

  /// damage
  CohesiveInternalField<Real> damage;

  /// pointer to the solid mechanics model for cohesive elements
  SolidMechanicsModelCohesive * model;

  /// critical stress
  RandomInternalField<Real, FacetInternalField> sigma_c;

  /// critical displacement
  Real delta_c;

  /// array to temporarily store the normals
  Array<Real> normal;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_cohesive_inline_impl.cc"

__END_AKANTU__

#include "cohesive_internal_field_tmpl.hh"

#endif /* __AKANTU_MATERIAL_COHESIVE_HH__ */
