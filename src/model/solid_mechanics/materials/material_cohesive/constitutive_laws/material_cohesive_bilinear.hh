/**
 * @file   material_cohesive_bilinear.hh
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Thu Feb 16 14:14:34 2012
 *
 * @brief  Bilinear cohesive constitutive law
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

#ifndef __AKANTU_MATERIAL_COHESIVE_BILINEAR_HH__
#define __AKANTU_MATERIAL_COHESIVE_BILINEAR_HH__

/* -------------------------------------------------------------------------- */

#include "material_cohesive_linear.hh"
#include "aka_common.hh"

__BEGIN_AKANTU__

/**
 * Cohesive material bilinear
 *
 * parameters in the material files :
 *   - delta_0   : elastic limit displacement (default: 0)
 */
template<UInt spatial_dimension>
class MaterialCohesiveBilinear : public MaterialCohesiveLinear<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  MaterialCohesiveBilinear(SolidMechanicsModel & model, const ID & id = "");
  virtual ~MaterialCohesiveBilinear();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize the material computed parameter
  virtual void initMaterial();

  /// update delta_max values with delta_0
  virtual void updateDeltaMax(GhostType ghost_type);

  /// resize vectors for new cohesive elements
  virtual void resizeCohesiveVectors();

protected:

  void computeTangentStiffness(__attribute__((unused))	const ElementType & el_type,
			       __attribute__((unused)) Vector<Real> & tangent_matrix,
			       __attribute__((unused)) GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
 // void computeTangentStiffness(__attribute__((unused)) Vector<Real> & tangent_matrix,
 // 				       __attribute__((unused)) const Vector<Real> & normal,
 // 			         	__attribute__((unused))	const ElementType & el_type,
 //  				       __attribute__((unused)) GhostType ghost_type = _not_ghost) {
 //    AKANTU_DEBUG_TO_IMPLEMENT();
 //  }



  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// elastic limit displacement
  Real delta_0;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "material_cohesive_elastic_inline_impl.cc"


__END_AKANTU__

#endif /* __AKANTU_MATERIAL_COHESIVE_ELASTIC_HH__ */
