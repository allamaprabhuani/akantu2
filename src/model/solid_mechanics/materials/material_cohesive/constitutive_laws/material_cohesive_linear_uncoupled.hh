/**
 * @file   material_cohesive_linear.hh
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Thu Jan 14 2016
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * random stress definition for extrinsic type
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

#include "material_cohesive_linear.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_COHESIVE_LINEAR_UNCOUPLED_HH__
#define __AKANTU_MATERIAL_COHESIVE_LINEAR_UNCOUPLED_HH__

/* -------------------------------------------------------------------------- */

namespace akantu {

/**
 * Cohesive material linear with two different laws for mode I and
 * mode II, for extrinsic case
 *
 * parameters in the material files :
 *  - roughness : define the interaction between mode I and mode II (default: 0)
 */
template<UInt spatial_dimension>
class MaterialCohesiveLinearUncoupled : public MaterialCohesiveLinear<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  //  typedef MaterialCohesiveLinear<spatial_dimension> MaterialParent;
public:

  MaterialCohesiveLinearUncoupled(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize the material parameters
  virtual void initMaterial();

protected:

  /// constitutive law
  virtual void computeTraction(const Array<Real> & normal,
			       ElementType el_type,
			       GhostType ghost_type = _not_ghost);

  /// check delta_max for cohesive elements in case of no convergence
  /// in the solveStep (only for extrinsic-implicit)
  virtual void checkDeltaMax(GhostType ghost_type = _not_ghost);

  /// compute tangent stiffness matrix
  virtual void computeTangentTraction(const ElementType & el_type,
				      Array<Real> & tangent_matrix,
				      const Array<Real> & normal,
				      GhostType ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// parameter to tune the interaction between mode II and mode I
  Real R;

  /// maximum normal displacement
  CohesiveInternalField<Real> delta_n_max;

  /// maximum tangential displacement
  CohesiveInternalField<Real> delta_t_max;

  /// damage associated to normal tractions 
  CohesiveInternalField<Real> damage_n;

  /// damage associated to shear tractions 
  CohesiveInternalField<Real> damage_t;

};

} // akantu

#endif /* __AKANTU_MATERIAL_COHESIVE_LINEAR_UNCOUPLED_HH__ */
