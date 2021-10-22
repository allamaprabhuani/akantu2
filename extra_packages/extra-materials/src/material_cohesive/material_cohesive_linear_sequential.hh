/**
 * @file   material_cohesive_linear.hh
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * random stress definition for extrinsic type
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* ------------------------------------------------------------------ */
#include "material_cohesive_linear.hh"
/* ------------------------------------------------------------------ */

#ifndef __AKANTU_MATERIAL_COHESIVE_LINEAR_SEQUENTIAL_HH__
#define __AKANTU_MATERIAL_COHESIVE_LINEAR_SEQUENTIAL_HH__

namespace akantu {

/**
 * Cohesive material linear damage for extrinsic case
 *
 * parameters in the material files :
 *   - sigma_c   : critical stress sigma_c  (default: 0)
 *   - beta      : weighting parameter for sliding and normal opening (default:
 * 0)
 *   - G_cI      : fracture energy for mode I (default: 0)
 *   - G_cII     : fracture energy for mode II (default: 0)
 *   - penalty   : stiffness in compression to prevent penetration
 */
template <UInt spatial_dimension>
class MaterialCohesiveLinearSequential
    : public MaterialCohesiveLinear<spatial_dimension> {
  /* ---------------------------------------------------------------- */
  /* Constructors/Destructors                                         */
  /* ---------------------------------------------------------------- */
public:
  MaterialCohesiveLinearSequential(SolidMechanicsModel & model,
                                   const ID & id = "");

  /* ---------------------------------------------------------------- */
  /* Methods                                                          */
  /* ---------------------------------------------------------------- */
public:
  /// initialise material
  void initMaterial() override;

  /// check stress for cohesive elements' insertion
  void checkInsertion(bool check_only = false) override;

  /// updates delta_max and damage (SLA functionality)
  UInt updateDeltaMax(ElementType el_type, GhostType ghost_type);

  /// identify a facet with the highest tensile stress complying with
  /// topological criterion
  std::tuple<UInt, Real, UInt>
  findCriticalFacet(const ElementType & type_facet);

  /// determine the crack contour_subfacet_coh_el, surface_subfacets_crac_nb,
  /// surface nodes and contour_nodes within this specific material
  std::tuple<std::map<Element, Element>, std::map<Element, UInt>,
             std::set<UInt>, std::set<UInt>>
  determineCrackSurface();

  /// searches for the "stressed" facets connected to contour segments
  std::map<UInt, UInt> findCriticalFacetsOnContour(
      const std::map<Element, Element> & contour_subfacets_coh_el,
      const std::map<Element, UInt> & surface_subfacets_crack_nb,
      const std::set<UInt> & contour_nodes,
      const std::set<UInt> & surface_nodes);

  /// update effective stresses, scalar and normal tractions on
  /// all facets of material
  void computeEffectiveStresses();

  /// insert a row of cohesives according to provided numbers
  void insertCohesiveElements(std::map<UInt, UInt> & facet_nbs_crack_nbs,
                              ElementType facet_type, bool check_only);

  /// insert facets having 2 or more cohesive neighbors
  std::map<UInt, UInt>
  findHolesOnContour(std::map<Element, Element> & contour_subfacets_coh_el,
                     const std::map<Element, UInt> & surface_subfacets_crack_nb,
                     const std::set<UInt> & surface_nodes);

  std::tuple<Real, Element, UInt>
  computeMaxDeltaMaxExcess(ElementType el_type,
                           GhostType ghost_type = _not_ghost);

protected:
  void computeTraction(const Array<Real> & normal, ElementType el_type,
                       GhostType ghost_type = _not_ghost) override;

  void computeTangentTraction(ElementType el_type, Array<Real> & tangent_matrix,
                              const Array<Real> & normal,
                              GhostType ghost_type) override;

  /// compute the traction based on a fixed stiffness (SLA style)
  inline void computeTractionOnQuad(
      Vector<Real> & traction, Vector<Real> & opening,
      const Vector<Real> & normal, Real & delta_max, const Real & delta_c,
      const Vector<Real> & insertion_stress, const Real & sigma_c,
      Vector<Real> & normal_opening, Vector<Real> & tangential_opening,
      Real & normal_opening_norm, Real & tangential_opening_norm, Real & damage,
      bool & penetration, Vector<Real> & contact_traction,
      Vector<Real> & contact_opening);

  /// compute traction with penalty = tensile stiffness
  inline void computeSimpleTractionOnQuad(
      Vector<Real> & traction, Vector<Real> & opening,
      const Vector<Real> & normal, Real & delta_max, const Real & delta_c,
      const Real & sigma_c, Vector<Real> & normal_opening,
      Vector<Real> & tangential_opening, Real & normal_opening_norm,
      Real & tangential_opening_norm, Real & damage, bool & penetration,
      Vector<Real> & contact_traction, Vector<Real> & contact_opening);

  /// compute the stiffness dependent only on previous delta max (SLA)
  inline void computeTangentTractionOnQuad(
      Matrix<Real> & tangent, Real & delta_max, const Real & delta_c,
      const Real & sigma_c, const Vector<Real> & normal,
      const Real & normal_opening_norm, const Real & tangential_opening_norm,
      const Real & damage);

  /// compute the stiffness dependent only on previous delta max both for
  /// tension and compression (SLA)
  inline void
  computeSecantTractionOnQuad(Matrix<Real> & tangent, Real & delta_max,
                              const Real & delta_c, const Real & sigma_c,
                              const Vector<Real> & normal,

                              const Real & damage, const Real & prev_damage);

  inline bool updateDeltaMaxOnQuad(const Real & normal_opening_norm,
                                   const Real & tangential_opening_norm,
                                   Real & damage, Real & delta_max,
                                   const Real & delta_c);

  inline Real computeDeltaMaxExcessOnQuad(const Real & normal_opening_norm,
                                          const Real & tangential_opening_norm,
                                          const Real & damage,
                                          const Real & delta_max,
                                          const Real & delta_c,
                                          bool & penetration);

  bool hasStiffnessMatrixChanged() override {
    UInt nb_element = 0;
    for (auto gt : ghost_types) {
      for (auto type : this->element_filter.elementTypes(spatial_dimension, gt,
                                                         _ek_cohesive)) {
        auto && elem_filter = this->element_filter(type, gt);
        nb_element += elem_filter.size();
      }
    }
    auto && comm = akantu::Communicator::getWorldCommunicator();
    comm.allReduce(nb_element, SynchronizerOperation::_sum);
    if (nb_element == 0) {
      return false;
    } else {
      comm.allReduce(update_stiffness, SynchronizerOperation::_lor);
      return (update_stiffness);
    }
  }

  /* ---------------------------------------------------------------- */
  /* Accessors                                                        */
  /* ---------------------------------------------------------------- */
public:
  /// get the effective stresses
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(EffectiveStress, effective_stresses,
                                         Real);

  AKANTU_SET_MACRO(UpdateStiffness, update_stiffness, bool);
  /* ---------------------------------------------------------------- */
  /* Class Members                                                    */
  /* ---------------------------------------------------------------- */
protected:
  /// scalar tractions = combination of normal and tangential norms
  FacetInternalField<Real> scalar_tractions;

  /// traction acting on a facet plane
  FacetInternalField<Real> normal_tractions;

  /// normal stresses normalized by the stress limit
  FacetInternalField<Real> effective_stresses;

  /// internal variable to indicate to solver if stiffness has to be reassembled
  bool update_stiffness{true};

  /// defines deviation from the exact value of the most stressed element to be
  /// damaged by formula delta < (1 + delta_deviation) * delta_max
  Real delta_deviation{0};

  /// number of maximum reductions within SLA
  UInt reductions{10};

  /// flag to update stiffness of an element if damaged
  bool update_stiffness_on_damage{true};
};

/* ------------------------------------------------------------------ */
/* inline functions                                                   */
/* ------------------------------------------------------------------ */

} // namespace akantu

#include "material_cohesive_linear_sequential_inline_impl.hh"

#endif /* __AKANTU_MATERIAL_COHESIVE_LINEAR_SEQUENTIAL_HH__ */
