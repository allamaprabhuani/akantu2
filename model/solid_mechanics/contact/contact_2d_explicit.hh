/**
 * @file   contact_2d_explicit.hh
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @date   Tue Oct 12 09:24:20 2010
 *
 * @brief  Interface for 2d explicit contact class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_2D_EXPLICIT_HH__
#define __AKANTU_CONTACT_2D_EXPLICIT_HH__

/* -------------------------------------------------------------------------- */

#include "aka_common.hh"
#include "contact.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */

class Contact2dExplicit : public Contact{
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Contact2dExplicit(const SolidMechanicsModel & model,
		    const ContactType & type,
		    const ContactID & id = "contact",
		    const MemoryID & memory_id = 0);

  virtual ~Contact2dExplicit();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void solveContact();

private:

  /// Project back nodes on penetrated segments
  void projectNodesOnSegements(PenetrationList & pen_list);

  /// Decompose velocities prior impact to compute normal velocities
  void computeNormalVelocities(PenetrationList & pen_list, Vector<Real> & gap_der, Vector<Real> & vel_norm);
  
  /// Decompose velocities prior impact to compute friction velocities
  void computeFrictionVelocities(PenetrationList & pen_list, Vector<Real> & gap_der, Vector<Real> & vel_norm, Vector<Real> & vel_fric);

  /// Update velocities adding normal and friction components
  void updatePostImpactVelocities(PenetrationList & pen_list,Vector<Real> & vel_norm, Vector<Real> & vel_fric);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// Set coefficient of restitution (default value = 0)
  AKANTU_SET_MACRO(CoefficientOfRestitution, coefficient_of_restitution, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// Coefficient of restitution used to compute post impact velocities
  Real coefficient_of_restitution;

};

__END_AKANTU__

#endif /* __AKANTU_CONTACT_2D_EXPLICIT_HH__ */
