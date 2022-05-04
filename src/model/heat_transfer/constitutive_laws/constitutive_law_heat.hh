/**
 * @file   constituive_law_heat.hh
 *
 * @author Mohit Pundir <mohit.pundir@ethz.ch>
 *
 * @date creation: Wed May 4 2022
 * @date last modification: Wed May 4 2022
 *
 * @brief  Heat transfer law for poisson model
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "constitutive_law.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONSTITUTIVE_LAW_HEAT_HH__
#define __AKANTU_CONSTITUTIVE_LAW_HEAT_HH__

namespace akantu {
class ConstitutiveLawHeat : public ConstitutiveLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ConstitutiveLawHeat(PoissonModel & model, const ID & id = "");

  ConstitutiveLawHeat(PoissonModel & model, UInt spatial_dimension,
		      const Mesh & mesh, FEEngine & fe_engine, const ID & id = "");

  
  ~ConstitutiveLawHeat() override = default;

protected:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial() override;

  /// constitutive law for all element of a type
  void computeFlux(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(ElementType el_type, Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;

  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// the density
  Real density;

  /// capacity
  Real capacity;

  // conductivity
  Real conductivity;

  std::unordered_map<GhostType, bool> initial_conductivity{{_not_ghost, true},
                                                           {_ghost, true}};

  std::unordered_map<GhostType, UInt> conductivity_release{{_not_ghost, 0},
                                                           {_ghost, 0}};
  
} // namespace akantu

#endif
