/**
 * @file   material_cohesive_linear.hh
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Mon Feb 20 12:00:34 2012
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading
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

#ifndef __AKANTU_MATERIAL_COHESIVE_LINEAR_HH__
#define __AKANTU_MATERIAL_COHESIVE_LINEAR_HH__

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/**
 * Cohesive material linear damage
 *
 * parameters in the material files :
 *   - sigma_c   : critical stress sigma_c  (default: 0)
 *   - beta      : weighting parameter for sliding and normal opening (default: 0)
 *   - G_cI      : fracture energy for mode I (default: 0)
 *   - G_cII     : fracture energy for mode II (default: 0)
 */

class MaterialCohesiveLinear : public MaterialCohesive {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  MaterialCohesiveLinear(Model & model, const ID & id = "");
  virtual ~MaterialCohesiveLinear();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// set patameters
  virtual bool setParam(const std::string & key, const std::string & value,
			const ID & id);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// initialize the material computed parameter
  virtual void initMaterial();

  /// resize vectors for new cohesive elements
  virtual void resizeCohesiveVectors();

protected:

  /// constitutive law
  void computeTraction(const Vector<Real> & normal,
		       ElementType el_type,
		       GhostType ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// critical stress
  Real sigma_c;

  /// beta parameter
  Real beta;

  /// mode I fracture energy
  Real G_cI;

  /// mode II fracture energy
  Real G_cII;

  /// kappa parameter
  Real kappa;

  /// maximum displacement
  ByElementTypeReal delta_max;

  /// critical displacement
  Real delta_c;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "material_cohesive_linear_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MaterialCohesiveLinear & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_MATERIAL_COHESIVE_LINEAR_HH__ */
