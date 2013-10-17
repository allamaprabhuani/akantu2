/**
 * @file   material_thermal.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date   Th Oct 17 11:48:20 2013
 *
 * @brief  Material isotropic thermo-elastic
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
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_THERMAL_HH__
#define __AKANTU_MATERIAL_THERMAL_HH__

__BEGIN_AKANTU__
template<UInt spatial_dimension>
class MaterialThermal : public virtual Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialThermal(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialThermal() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();
  void updateInternalParameters();

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Thermal expansion coefficient
  Real alpha;

  /// Temperature field
  ByElementTypeReal delta_t;

  /// Uniform temperature field (for test puposes)
  Real delta_t_init;
};

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_THERMAL_HH__ */
