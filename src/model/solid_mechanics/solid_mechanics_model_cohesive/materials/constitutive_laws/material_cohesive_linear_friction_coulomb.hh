/**
 * @file   material_cohesive_linear_friction_coulomb.hh
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Mathias Lebihain <mathias.lebihain@epfl.ch>
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

/* -------------------------------------------------------------------------- */

#include "material_cohesive_linear_friction.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_COHESIVE_LINEAR_FRICTION_COULOMB_HH__
#define __AKANTU_MATERIAL_COHESIVE_LINEAR_FRICTION_COULOMB_HH__

/* -------------------------------------------------------------------------- */

namespace akantu {

/**
 * Cohesive material linear with coulomb friction force
 */
template <UInt spatial_dimension>
class MaterialCohesiveLinearFrictionCoulomb
    : public MaterialCohesiveLinearFriction<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  using MaterialParent = MaterialCohesiveLinearFriction<spatial_dimension>;

public:
  MaterialCohesiveLinearFrictionCoulomb(SolidMechanicsModel & model,
                                        const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// constitutive law
  void computeTraction(const Array<Real> & normal, ElementType el_type,
                       GhostType ghost_type = _not_ghost) override;

  // /// compute tangent stiffness matrix
  // void computeTangentTraction(const ElementType & el_type,
  //                             Array<Real> & tangent_matrix,
  //                             const Array<Real> & normal,
  //                             GhostType ghost_type) override;
};

} // namespace akantu

#endif /* __AKANTU_MATERIAL_COHESIVE_LINEAR_FRICTION_COULOMB_HH__ */
