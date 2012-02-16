/**
 * @file   velocity_dependent_fric_coef.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Jun 20 15:19:49 2011
 *
 * @brief  implementation of velocity dependence for friction coefficient
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

#ifndef __AKANTU_VELOCITY_DEPENDENT_FRIC_COEF_HH__
#define __AKANTU_VELOCITY_DEPENDENT_FRIC_COEF_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "contact_rigid.hh"
#include "friction_coefficient.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class VelocityDependentFricCoef : public virtual FrictionCoefficient {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  VelocityDependentFricCoef(ContactRigid & contact,
			    const Surface & master_surface);

  virtual ~VelocityDependentFricCoef();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// compute the relative sliding velocity
  virtual void initializeComputeFricCoef();

  /// implementation of the friction coefficient formula
  virtual Real computeFricCoef(UInt impactor_node_index) = 0;

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;
  
private:
  // computes the tangential velocity of the master element
  void computeTangentialMasterVelocity(UInt impactor_index, 
				       ContactRigid::ImpactorInformationPerMaster * impactor_info, 
				       Real * tangential_master_velocity);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// spatial dimension of contact
  UInt spatial_dimension;

  /// relative sliding velocities for each active impactor node
  Vector<Real> * relative_sliding_velocities;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "velocity_dependent_fric_coef_inline_impl.cc"

/*
/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const VelocityDependentFricCoef & _this)
{
  _this.printself(stream);
  return stream;
}
*/


__END_AKANTU__

#endif /* __AKANTU_VELOCITY_DEPENDENT_FRIC_COEF_HH__ */
