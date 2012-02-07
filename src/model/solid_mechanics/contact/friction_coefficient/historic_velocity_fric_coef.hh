/**
 * @file   historic_velocity_fric_coef.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Nov  1 10:48:22 2011
 *
 * @brief  friction coefficient that depends on the historic of velocity
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

#ifndef __AKANTU_HISTORIC_VELOCITY_FRIC_COEF_HH__
#define __AKANTU_HISTORIC_VELOCITY_FRIC_COEF_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_circular_vector.hh"
#include "contact_rigid.hh"
#include "friction_coefficient.hh"

__BEGIN_AKANTU__

class HistoricVelocityFricCoef : public virtual FrictionCoefficient {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /// \warning DO NOT USE THIS CLASS !!!

  HistoricVelocityFricCoef(ContactRigid & contact,
			   const Surface & master_surface,
			   const Real beta);

  // constructur when this class is not used
  HistoricVelocityFricCoef(ContactRigid & contact,
			   const Surface & master_surface);
  
  virtual ~HistoricVelocityFricCoef();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// compute the history of relative sliding velocity
  virtual void initializeComputeFricCoef();

  /// implementation of the friction coefficient formula
  virtual Real computeFricCoef(UInt impactor_node_index) = 0;
  //  virtual Real computeFricCoef(UInt impactor_node_index) { return 0.;}; // for testing
  
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
private:
  /// spatial dimension of contact
  UInt spatial_dimension;

protected:
  /// weight parameter
  Real beta;
  
  /// weights
  Vector<Real> * weights;

  /// history of sliding velocities
  CircularVector< Vector<Real> * > * historic_velocities;
  
  /// order of information
  Vector<UInt> nodes;

  /// if node is active in this time step
  Vector<bool> active;

  /// time since node came in contact
  Vector<Real> contact_time;

  /// keep information about stick status of active impactor node
  Vector<bool> * node_stick_status;
  
  //public: // for testing
  /// relative sliding velocities for each active impactor node
  Vector<Real> * generalized_sliding_velocities;
};


/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "historic_velocity_fric_coef_inline_impl.cc"

/*
/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const HistoricVelocityFricCoef & _this)
{
  _this.printself(stream);
  return stream;
}
*/

__END_AKANTU__


#endif /* __AKANTU_HISTORIC_VELOCITY_FRIC_COEF_HH__ */
