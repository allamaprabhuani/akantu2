/**
 * @file   constitutive_law_non_local.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Mon Sep 11 2017
 *
 * @brief  ConstitutiveLaw class that handle the non locality of a law for
 * example damage.
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

/* -------------------------------------------------------------------------- */
#include "constitutive_law.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_CONSTITUTIVE_LAW_NON_LOCAL_HH_
#define AKANTU_CONSTITUTIVE_LAW_NON_LOCAL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
class ConstitutiveLawNonLocalInterface {
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the constitutive_law the non local parts of the
  /// constitutive_law
  void initConstitutiveLawNonLocal() {
    this->registerNeighborhood();
    this->registerNonLocalVariables();
  };

  /// insert the quadrature points in the neighborhoods of the non-local manager
  virtual void insertIntegrationPointsInNeighborhoods(
      GhostType ghost_type,
      const ElementTypeMapReal & quadrature_points_coordinates) = 0;

  /// update the values in the non-local internal fields
  virtual void updateNonLocalInternals(ElementTypeMapReal & non_local_flattened,
                                       const ID & field_id,
                                       GhostType ghost_type,
                                       ElementKind kind) = 0;

protected:
  /// get the name of the neighborhood for this constitutive_law
  virtual ID getNeighborhoodName() = 0;

  /// register the neighborhoods for the constitutive_law
  virtual void registerNeighborhood() = 0;

  /// register the non local internal variable
  virtual void registerNonLocalVariables() = 0;

  virtual inline void onElementsAdded(const Array<Element> & /*unused*/,
                                      const NewElementsEvent & /*unused*/) {}
};
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <UInt dim, class ConstitutiveLawNonLocalInterface, class ConstitutiveLawParent>
class ConstitutiveLawNonLocal : public ConstitutiveLawNonLocalInterface,
                                public ConstitutiveLawParent {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  explicit ConstitutiveLawNonLocal(
      typename ConstitutiveLawParent::ConstitutiveLawsHandler & handler,
      const ID & id);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// insert the quadrature points in the neighborhoods of the non-local manager
  void insertIntegrationPointsInNeighborhoods(
      GhostType ghost_type,
      const ElementTypeMapReal & quadrature_points_coordinates) override;

  /// update the values in the non-local internal fields
  void updateNonLocalInternals(ElementTypeMapReal & non_local_flattened,
                               const ID & field_id, GhostType ghost_type,
                               ElementKind kind) override;

  /// register the neighborhoods for the constitutive_law
  void registerNeighborhood() override;

protected:
  /// get the name of the neighborhood for this constitutive_law
  ID getNeighborhoodName() override { return this->name; }
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
#include "constitutive_law_non_local_interface_tmpl.hh"

#endif /* AKANTU_CONSTITUTIVE_LAW_NON_LOCAL_HH_ */
