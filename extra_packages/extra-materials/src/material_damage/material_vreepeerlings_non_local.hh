/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material_damage_non_local.hh"
#include "material_vreepeerlings.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_VREEPEERLINGS_NON_LOCAL_HH_
#define AKANTU_MATERIAL_VREEPEERLINGS_NON_LOCAL_HH_

namespace akantu {

/**
 * Material VreePeerlings Non local
 *
 * parameters in the material files :
 */
template <Int spatial_dimension,
          template <UInt> class MatParent = MaterialElastic>
class MaterialVreePeerlingsNonLocal
    : public MaterialDamageNonLocal<
          spatial_dimension,
          MaterialVreePeerlings<spatial_dimension, MatParent>> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef MaterialVreePeerlings<spatial_dimension, MatParent> Parent;
  typedef MaterialDamageNonLocal<spatial_dimension, Parent>
      MaterialVreePeerlingsNonLocalParent;

  MaterialVreePeerlingsNonLocal(SolidMechanicsModel & model,
                                const ID & id = "");

  virtual ~MaterialVreePeerlingsNonLocal(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial();

  /// constitutive law for all element of a type
  // void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// constitutive law
  virtual void computeNonLocalStress(ElementType el_type,
                                     GhostType ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// non local version of equivalent strain
  InternalField<Real> equi_strain_non_local;

  /// non local version of equivalent strain rate
  InternalField<Real> equi_strain_rate_non_local;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_vreepeerlings_non_local_inline_impl.hh"

} // namespace akantu

#endif /* AKANTU_MATERIAL_VREEPEERLINGS_NON_LOCAL_HH_ */
