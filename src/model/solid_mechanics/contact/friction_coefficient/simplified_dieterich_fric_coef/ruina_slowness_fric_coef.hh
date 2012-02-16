/**
 * @file   ruina_slowness_fric_coef.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Mar  7 16:08:22 2011
 *
 * @brief  the friction law based on the slowness law of the theta state variable
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

#ifndef __AKANTU_RUINA_SLOWNESS_FRIC_COEF_HH__
#define __AKANTU_RUINA_SLOWNESS_FRIC_COEF_HH__

/* -------------------------------------------------------------------------- */
#include "contact_rigid.hh"
#include "simplified_dieterich_fric_coef.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

template<bool compute_analytic_solution = true>
class RuinaSlownessFricCoef : public SimplifiedDieterichFricCoef {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  RuinaSlownessFricCoef(ContactRigid & contact,
			const Surface & master_surface);

  virtual ~RuinaSlownessFricCoef();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void setParam(const std::string & key, const std::string & value);

  /// compute the relative sliding velocity and the theta state variable
  virtual void initializeComputeFricCoef();

  /// add an impactor surface to this master surface
  virtual void addImpactorSurface(const Surface & impactor_surface);

  /// remove an impactor surface of this master surface
  virtual void removeImpactorSurface(const Surface & impactor_surface);

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  /// compute the state variable by the analytical solution
  inline Real computeAnalyticTheta(Real previous_theta, Real sliding_speed, Real delta_t);

  /// compute the state variable by the implicit solution
  inline Real computeImplicitTheta(Real previous_theta, Real sliding_speed, Real delta_t);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// compute the analytic solution of theta or the implicit solution
  //  bool compute_analytic_solution;

  /// the characteristic length d_0
  Real d_zero;

  /// assigns the index to an impactor node
  std::map<UInt, UInt> node_to_index;

  /// the active impactor nodes from the precedent time step
  Vector<bool> * are_active_impactor_nodes;
  Vector<bool> * were_active_impactor_nodes;

  /// the internal theta state variable
  Vector<Real> * previous_theta_state_variables;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "ruina_slowness_fric_coef_inline_impl.cc"
#endif
/*
/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const RuinaSlownessFricCoef & _this)
{
  _this.printself(stream);
  return stream;
}
*/

__END_AKANTU__

#endif /* __AKANTU_RUINA_SLOWNESS_FRIC_COEF_HH__ */
