/**
 * @file   material_non_local.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Dec 08 2015
 *
 * @brief  Material class that handle the non locality of a law for example
 * damage.
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
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_NON_LOCAL_HH__
#define __AKANTU_MATERIAL_NON_LOCAL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
class MaterialNonLocalInterface {
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material computed parameter
  virtual void initMaterial() { this->registerNeighborhood(); }

  /// insert the quadrature points in the neighborhoods of the non-local manager
  virtual void insertIntegrationPointsInNeighborhoods(
      const GhostType & ghost_type,
      const ElementTypeMapReal & quadrature_points_coordinates) = 0;

  /// update the values in the non-local internal fields
  virtual void updateNonLocalInternals(ElementTypeMapReal & non_local_flattened,
                                       const ID & field_id,
                                       const GhostType & ghost_type,
                                       const ElementKind & kind) = 0;
  /// constitutive law
  virtual void computeNonLocalStresses(GhostType ghost_type = _not_ghost) = 0;

  /// register the neighborhoods for the material
  virtual void registerNeighborhood() = 0;

protected:
  virtual inline void onElementsAdded(const Array<Element> &,
                                      const NewElementsEvent &) {}
};
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <UInt dim, class LocalParent>
class MaterialNonLocal : public MaterialNonLocalInterface, public LocalParent {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  explicit MaterialNonLocal(SolidMechanicsModel & model, const ID & id);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material computed parameter
  void initMaterial() override;

  /// insert the quadrature points in the neighborhoods of the non-local manager
  void insertIntegrationPointsInNeighborhoods(
      const GhostType & ghost_type,
      const ElementTypeMapReal & quadrature_points_coordinates) override;

  /// update the values in the non-local internal fields
  void updateNonLocalInternals(ElementTypeMapReal & non_local_flattened,
                               const ID & field_id,
                               const GhostType & ghost_type,
                               const ElementKind & kind) override;

  /// register the neighborhoods for the material
  void registerNeighborhood() override;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
#include "material_non_local_tmpl.hh"

#endif /* __AKANTU_MATERIAL_NON_LOCAL_HH__ */
