/**
 * @file   integrator.hh
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Thu Feb 10 11:09:12 2011
 *
 * @brief  interface for integrator classes
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

#ifndef __AKANTU_INTEGRATOR_HH__
#define __AKANTU_INTEGRATOR_HH__

/* -------------------------------------------------------------------------- */
#include "aka_memory.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

class Integrator : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Integrator(Mesh & m, IntegratorID myid="integrator"){
    mesh = &m;
    id = myid;
  };
  virtual ~Integrator(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  template <ElementType type>
  inline void precomputeJacobiansOnQuadraturePoints(__attribute__ ((unused)) const UInt dimension,
						    __attribute__ ((unused)) GhostType ghost_type) { }

  void integrateOnElement(__attribute__ ((unused)) const Vector<Real> & f,
			  __attribute__ ((unused)) Real * intf,
			  __attribute__ ((unused)) UInt nb_degre_of_freedom,
			  __attribute__ ((unused)) const Element & elem,
			  __attribute__ ((unused)) GhostType ghost_type) const {};

  /// function to print the contain of the class
  //  virtual void printself(std::ostream & stream, int indent = 0) const{};

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */


protected:

  Mesh * mesh;
  IntegratorID id;


  /// jacobians for all elements
  ByElementTypeReal jacobians;

  /// jacobians for all elements
  ByElementTypeReal ghost_jacobians;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "integrator_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const Integrator & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__

#endif /* __AKANTU_INTEGRATOR_HH__ */
