/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "material.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include "cohesive_internal_field.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_COHESIVE_HH_
#define AKANTU_MATERIAL_COHESIVE_HH_

/* -------------------------------------------------------------------------- */
namespace akantu {
class SolidMechanicsModelCohesive;
}

namespace akantu {

class MaterialCohesive : public Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using MyFEEngineCohesiveType =
      FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_cohesive>;

public:
  MaterialCohesive(SolidMechanicsModel & model, const ID & id = "");
  ~MaterialCohesive() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material computed parameter
  void initMaterial() override;

  /// compute tractions (including normals and openings)
  void computeTraction(GhostType ghost_type = _not_ghost);

  /// assemble residual
  void assembleInternalForces(GhostType ghost_type = _not_ghost) override;

  /// check stress for cohesive elements' insertion, by default it
  /// also updates the cohesive elements' data
  virtual void checkInsertion(bool /*check_only*/ = false) {
    AKANTU_TO_IMPLEMENT();
  }

  /// interpolate   stress  on   given   positions  for   each  element   (empty
  /// implemantation to avoid the generic call to be done on cohesive elements)
  virtual void interpolateStress(const ElementType /*type*/,
                                 Array<Real> & /*result*/) {}

  /// compute the stresses
  void computeAllStresses(GhostType /*ghost_type*/ = _not_ghost) final{};

  // add the facet to be handled by the material
  Idx addFacet(const Element & element);

  template <Int dim>
  inline decltype(auto) getArguments(ElementType element_type,
                                     GhostType ghost_type) {
    return zip(
        "normal"_n = make_view<dim>(this->normals(element_type, ghost_type)),
        "opening"_n = make_view<dim>(this->opening(element_type, ghost_type)),
        "traction"_n =
            make_view<dim>(this->tractions(element_type, ghost_type)),
        "contact_opening"_n =
            make_view<dim>(this->contact_opening(element_type, ghost_type)),
        "contact_traction"_n =
            make_view<dim>(this->contact_tractions(element_type, ghost_type)),
        "previous_opening"_n =
            make_view<dim>(this->opening.previous(element_type, ghost_type)),
        "previous_traction"_n =
            make_view<dim>(this->tractions.previous(element_type, ghost_type)),
        "delta_max"_n = this->delta_max(element_type, ghost_type),
        "damage"_n = this->damage(element_type, ghost_type));
  }

protected:
  virtual void computeTangentTraction(ElementType /*el_type*/,
                                      Array<Real> & /*tangent_matrix*/,
                                      GhostType /*ghost_type*/ = _not_ghost) {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute the normal
  void computeNormal(const Array<Real> & position, Array<Real> & normal,
                     ElementType type, GhostType ghost_type);

  /// compute the opening
  void computeOpening(const Array<Real> & displacement, Array<Real> & opening,
                      ElementType type, GhostType ghost_type);

  template <ElementType type>
  void computeNormal(const Array<Real> & position, Array<Real> & normal,
                     GhostType ghost_type);

  /// assemble stiffness
  void assembleStiffnessMatrix(GhostType ghost_type) override;

  /// constitutive law
  virtual void computeTraction(ElementType el_type,
                               GhostType ghost_type = _not_ghost) = 0;

  /// parallelism functions
  inline Int getNbData(const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

protected:
  void updateEnergies(ElementType el_type) override;

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
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(FacetFilter, facet_filter, Idx);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(FacetFilter, facet_filter, Idx);
  AKANTU_GET_MACRO_AUTO(FacetFilter, facet_filter);
  // AKANTU_GET_MACRO(ElementFilter, element_filter, const
  // ElementTypeMapArray<UInt> &);

  /// compute reversible energy
  auto getReversibleEnergy() -> Real;

  /// compute dissipated energy
  auto getDissipatedEnergy() -> Real;

  /// compute contact energy
  auto getContactEnergy() -> Real;

  /// get energy
  auto getEnergy(const std::string & type) -> Real override;

  /// return the energy (identified by id) for the provided element
  // Real getEnergy(const std::string & energy_id, const Element & element
  //                ) override {
  //   return Material::getEnergy(energy_id, element);
  // }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// list of facets assigned to this material
  ElementTypeMapArray<Idx> facet_filter;

  /// Link to the cohesive fem object in the model
  FEEngine & fem_cohesive;

private:
  /// reversible energy by quadrature point
  CohesiveInternalField<Real> reversible_energy;

  /// total energy by quadrature point
  CohesiveInternalField<Real> total_energy;

protected:
  /// opening in all elements and quadrature points
  CohesiveInternalField<Real> opening;

  /// traction in all elements and quadrature points
  CohesiveInternalField<Real> tractions;

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
  CohesiveInternalField<Real> normals;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "cohesive_internal_field_tmpl.hh"
#include "material_cohesive_inline_impl.hh"

#endif /* AKANTU_MATERIAL_COHESIVE_HH_ */
