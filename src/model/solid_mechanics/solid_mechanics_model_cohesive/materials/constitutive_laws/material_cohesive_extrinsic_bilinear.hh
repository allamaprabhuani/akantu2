/**
 * @file   material_cohesive_linear_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Shenghan Zhang <shenghan.zhang@epfl.ch>
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Dec 1 2016
 * @date last modification: Fri Mar 5 2021
 *
 * @brief  Bilinear Cohesive cohesive law for extrinsic elements
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_cohesive_linear.hh"

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_COHESIVE_EXTRINSIC_BILINEAR_HH__
#define __AKANTU_MATERIAL_COHESIVE_EXTRINSIC_BILINEAR_HH__

/* -------------------------------------------------------------------------- */

namespace akantu {

template<UInt spatial_dimension>
class MaterialCohesiveExtrinsicBilinear
  : public MaterialCohesiveLinear<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialCohesiveExtrinsicBilinear(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material computed parameter 
  void initMaterial() override;
  

protected:
  /// constitutive law
  void computeTraction(const Array<Real> & normal, ElementType el_type,
		       GhostType ghost_type = _not_ghost) override;
  
  /// compute the traction for a given quadrature point
  inline void computeTractionOnQuad(
    Vector<Real> & traction, Vector<Real> & opening,
    const Vector<Real> & normal, Real & delta_max, const Real & delta_c,
    const Real & delta_h, const Vector<Real> & insertion_stress,
    const Real & sigma_c, Vector<Real> & normal_opening,
    Vector<Real> & tangential_opening,
    Real & normal_opening_norm, Real & tangential_opening_norm, Real & damage,
    bool & penetration, Vector<Real> & contact_traction,
    Vector<Real> & contact_opening);


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  
  Real h;

  Real delta_h;
  
  CohesiveInternalField<Real> delta_h_eff;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */


}  // namespace akantu

#include "material_cohesive_extrinsic_bilinear_inline_impl.hh"


#endif /* __AKANTU_MATERIAL_COHESIVE_EXTRINSIC_BILINEAR_HH__ */
