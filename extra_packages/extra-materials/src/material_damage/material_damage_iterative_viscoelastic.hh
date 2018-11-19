/**
 * @file
 * @author
 * @date
 *
 * @brief
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
#include "material_iterative_stiffness_reduction.hh"
#include "material_viscoelastic_maxwell.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_DAMAGE_ITERATIVE_VISCOELASTIC_HH__
#define __AKANTU_MATERIAL_DAMAGE_ITERATIVE_VISCOELASTIC_HH__

namespace akantu {

/**
 * Material damage iterative viscoelastic
 */

/* -------------------------------------------------------------------------- */
template<UInt dim>
class MaterialDamageIterativeViscoelastic : public MaterialIterativeStiffnessReduction<dim, MaterialViscoelasticMaxwell> {
public:
  MaterialDamageIterativeViscoelastic(SolidMechanicsModel & model,
                                      const ID & id = "") : MaterialIterativeStiffnessReduction<dim, MaterialViscoelasticMaxwell>(model, id) {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// compute the tangent stiffness matrix for an element type
  void computeTangentModuli(const ElementType & el_type,
                            Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override;

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */

public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// compute tangent moduli on a quadrature point
  void computeTangentModuliOnQuad(Matrix<Real> & tangent, Real & dam);

};

} // namespace akantu

#include "material_damage_iterative_viscoelastic_inline_impl.cc"

#endif /* __AKANTU_MATERIAL_DAMAGE_ITERATIVE_VISCOELASTIC_HH__ */
