/**
 * @file   simplified_dieterich_fric_coef.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Fri Mar  4 14:51:28 2011
 *
 * @brief  implementation of a friction coefficient based on the simpliefied dieterich law \mu = \mu_0 + A \ln\left( \frac{V}{V^*} + 1\right) + B \ln \left( \frac{\theta}{\theta^*} + 1 \right) (see Dieterich & Kilgore 1994)
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

#ifndef __AKANTU_SIMPLIFIED_DIETERICH_FRIC_COEF_HH__
#define __AKANTU_SIMPLIFIED_DIETERICH_FRIC_COEF_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "contact_rigid.hh"
#include "friction_coefficient.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class SimplifiedDieterichFricCoef : public FrictionCoefficient {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  SimplifiedDieterichFricCoef(ContactRigid & contact,
			      const Surface & master_surface);

  virtual ~SimplifiedDieterichFricCoef();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void setParam(const std::string & key, const std::string & value);

  /// compute the relative sliding velocity
  virtual void initializeComputeFricCoef();

  /// implementation of the friction coefficient formula
  __aka_inline__ Real computeFricCoef(UInt impactor_node_index);

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
  Real spatial_dimension;

  /// /mu_0, the reference friction coefficient
  Real mu_zero;
  
  /// A, the factor of the rate term
  Real a_factor;

  /// B, the factor of the state term
  Real b_factor;

  /// The velocity normalizing constant
  Real v_normalizer;

  /// The state normalizing constant
  Real theta_normalizer;

  /// relative sliding velocities for each active impactor node
  Vector<Real> * relative_sliding_velocities;

  /// state for each active impactor node
  Vector<Real> * theta_state_variables;
};


/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "simplified_dieterich_fric_coef_inline_impl.cc"
#endif

/*
/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const SimplifiedDieterichFricCoef & _this)
{
  _this.printself(stream);
  return stream;
}
*/


__END_AKANTU__

#endif /* __AKANTU_SIMPLIFIED_DIETERICH_FRIC_COEF_HH__ */
