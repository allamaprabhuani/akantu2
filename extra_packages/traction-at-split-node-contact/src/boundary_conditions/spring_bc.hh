/**
 * @file   spring_bc.hh
 *
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  spring boundary condition
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __AST_SPRING_BC_HH__
#define __AST_SPRING_BC_HH__

// simtools
#include "force_based_dirichlet.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
class SpringBC : public ForceBasedDirichlet {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SpringBC(SolidMechanicsModel & model, 
	   BC::Axis ax, 
	   Real stiffness,
	   Real mass = 0.) :
    ForceBasedDirichlet(model, ax, 0., mass), 
    stiffness(stiffness), 
    elongation(0.)
  {}

  virtual ~SpringBC() {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual Real update() {
    AKANTU_DEBUG_IN();

    this->target_force = - this->stiffness * this->elongation;
    Real reaction = ForceBasedDirichlet::update();
    this->elongation += this->value;

    AKANTU_DEBUG_OUT();
    return reaction;
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Elongation, elongation, Real);

  inline void setToEquilibrium() {
    AKANTU_DEBUG_IN();
  
    this->updateTotalResidual();

    this->target_force = - this->total_residual;
    this->elongation = - this->target_force / this->stiffness;
  
    AKANTU_DEBUG_OUT();
  }

  /// change elongation
  /// dx > 0 -> target_force < 0
  inline void incrementElongation(Real dx) {
    AKANTU_DEBUG_IN();

    this->elongation += dx;

    AKANTU_DEBUG_OUT();
  }
  //friend std::ostream& operator<<(std::ostream& out, const SpringBC & spring);
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  Real stiffness;
  Real elongation;
};










// class SpringBCRestricted : public SpringBC {
// public:
//   SpringBCRestricted(BC::Axis ax, Real target_force, BC::Axis surface_axis, Real min, Real max)
//     :SpringBC(ax, target_force), surface_axis(surface_axis), min(min), max(max) {}

//   virtual ~SpringBCRestricted() {}

// public:
//   inline void operator()(UInt node, Vector<bool> & flags, Vector<Real> & primal, const Vector<Real> & coord) const {
//     if(coord(surface_axis) > min && coord(surface_axis) < max) {
//       SpringBC::operator()(node, flags, primal, coord);
//     }
//   }
// private:
//   BC::Axis surface_axis;
//   Real min;
//   Real max;
// };

// std::ostream& operator<<(std::ostream& out, const SpringBC & spring) {
//   out << "Real total_residual: " << *spring.total_residual << std::endl;
//   out << "Real mass: " << spring.mass << std::endl;
//   out << "Real k: " << spring.k << std::endl;
//   out << "Real delta_x: " << spring.delta_x << std::endl;
//   out << "Real dt: " << spring.dt << std::endl;
//   out << "Real v: " << spring.v << std::endl;
//   out << "Real dx: " << spring.dx << std::endl;
//   out << "Real forcing_vel: " << spring.forcing_vel << std::endl;
//   return out;
// }


__END_AKANTU__

#endif /* __AST_SPRING_BC_HH__ */
