/**
 * @file   material_cohesive_exponential.hh
 *
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Mon Jul 09 14:13:56 2012
 *
 * @brief  Exponential irreversible cohesive law of mixed mode loading
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
#include "material_cohesive.hh"
#include "aka_common.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_COHESIVE_EXPONENTIAL_HH__
#define __AKANTU_MATERIAL_COHESIVE_EXPONENTIAL_HH__

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/**
 * Cohesive material Exponential damage
 *
 * parameters in the material files :
 *   - sigma_c   : critical stress sigma_c  (default: 0)
 *   - beta      : weighting parameter for sliding and normal opening (default: 0)
 *   - delta_c   : critical opening (default: 0)
 */
template<UInt spatial_dimension>
class MaterialCohesiveExponential : public MaterialCohesive {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialCohesiveExponential(SolidMechanicsModel & model, const ID & id = "");
  virtual ~MaterialCohesiveExponential();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize the material computed parameter
  virtual void initMaterial();

  /// resize vectors for new cohesive elements
  virtual void resizeCohesiveVectors();



protected:

  /// constitutive law
  void computeTraction(const Vector<Real> & normal,
		       ElementType el_type,
		       GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  void computeTangentTraction(const ElementType & el_type,
			      Vector<Real> & tangent_matrix,
			      const Vector<Real> & normal,
			      GhostType ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:


  /// beta parameter
  Real beta;

  // /// maximum displacement
  // ByElementTypeReal delta_max;

  /// critical displacement
  Real delta_c;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_cohesive_exponential_inline_impl.cc"


__END_AKANTU__

#endif /* __AKANTU_MATERIAL_COHESIVE_EXPONENTIAL_HH__ */
