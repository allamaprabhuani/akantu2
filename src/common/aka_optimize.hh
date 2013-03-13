/**
 * @file   aka_optimize.hh
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @date   Fri Apr 20 10:44:00 2012
 *
 * @brief  Objects that can be used to carry out optimization
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

#ifndef __AKANTU_OPTIMIZE_HH__
#define __AKANTU_OPTIMIZE_HH__

#include <iostream>


#include <nlopt.hpp>
#include <array/expr.hpp>

#include "aka_config.hh"
#include "aka_common.hh"
#include "aka_point.hh"
#include "solid_mechanics_model.hh"

//#define DEBUG_OPTIMIZE 1


__BEGIN_AKANTU__

using std::cout;
using std::endl;

typedef array::Array<1,Real> vector_type;
typedef array::Array<2,Real> matrix_type;

std::ostream& operator<<(std::ostream&, nlopt::result);



enum Optimizator_type { Min_t, Max_t };


class Optimizator : public nlopt::opt {
  
  typedef std::vector<Real> point_type;
  typedef nlopt::opt base_type;
  
  point_type& x_;
  Real min_;
  int &count_;

public:
    
  // functor parameter constructor
  template <class functor_type>
  Optimizator(point_type& x0,
              functor_type& fn,
              Optimizator_type t = Min_t,
              nlopt::algorithm alg = nlopt::LD_SLSQP) : 
  base_type(alg, x0.size()), x_(x0), count_(fn.counter()) {
    
    if (t == Min_t)
      this->set_min_objective(functor_type::wrap, &fn);
    else
      this->set_max_objective(functor_type::wrap, &fn);
    
    this->set_xtol_rel(1e-4);
  }
  
  Real result() {
    
    optimize(x_, min_);
    
    cout<<"*** INFO *** Optimum value found at location";
    for (size_t i=0; i<x_.size(); ++i)
      cout<<" "<<x_[i];
    cout<<" after " <<count_<<" evaluations: "<<min_<<endl;
    
    return min_;
  }
  
};



template <ElementType>
class Distance_minimzer;

template <>
class Distance_minimzer<_segment_2> {
  
  static const UInt d = 2;
  static const UInt nb_nodes = 2;
  
  typedef Point<d> point_type;
  
  nlopt::opt opt_;         //!< Optimizator reference
  std::vector<Real> xi_;   //!< Master coordinate closest to point
  vector_type p_;          //!< Point to which the distance is minimized
  matrix_type XX_;         //!< Triangle coordinates
  UInt counter_;           //!< Optimization iteration counter
  Real fmin_;              //!< Minimum distance value
  
public:
  
  Distance_minimzer(const Real *r, const Element *el, SolidMechanicsModel &model)
  : opt_(nlopt::LD_SLSQP, d-1), xi_(d-1), p_(d), XX_(nb_nodes,d), counter_() {
    
    // set optimization parameters
    std::vector<Real> lb(1,-1), ub(1,1);
    opt_.set_lower_bounds(lb);
    opt_.set_upper_bounds(ub);
    opt_.set_min_objective(wrap, this);
    opt_.set_ftol_abs(1e-4);
    
    Mesh& mesh = model.getMesh();
    const Array<Real> &X = model.getCurrentPosition();
    const Array<UInt> &conn = mesh.getConnectivity(el->type);
    for (UInt i=0; i<nb_nodes; ++i) {
      XX_(0u,i) = X(conn(el->element,0),i);
      XX_(1u,i) = X(conn(el->element,1),i);
      p_(i) = r[i];
    }
    
    // compute start point
    start();
    
    // optimize
#ifdef DEBUG_OPTIMIZE
    nlopt::result result = opt_.optimize(xi_, fmin_);
    if (result > 0)
      cout<<"Optimium found in "<<counter_<<" iterations: "<<fmin_<<endl;
    cout<<"Point at master coordinate "<<xi_[0]<<": "<<point()<<endl;
    cout<<result<<endl;
#else
    opt_.optimize(xi_, fmin_);
#endif
  }
  
  vector_type operator()(const std::vector<Real> &xi)
  {
    vector_type N(nb_nodes);
    ElementClass<_segment_2>::computeShapes(&xi[0], 1, &N(0));
    return transpose(XX_)*N;
  }
  
  void start() {
    
    Real min = std::numeric_limits<Real>::infinity();
    Real xstart[3] = { -1., 0., 1. }; // check center and extremes of element
    int idx;
    for (int i=0; i<3; ++i) {
      xi_[0] = xstart[i];
      std::vector<Real> grad; // empty vector
      Real new_dist = (*this)(xi_, grad);
      if (new_dist < min) {
        min = new_dist;
        idx = i;
      }
    }
    xi_[0] = xstart[idx];
  }
  
  
  Real operator()(const std::vector<Real> &xi, std::vector<Real> &grad)
  {
    // increment function evaluation counter
    ++counter_;
    
    vector_type x = (*this)(xi);
    vector_type diff = x-p_;
    
    if (!grad.empty()) {
      
      // compute shape function derivatives
      matrix_type DN(nb_nodes, d-1);
      ElementClass<_segment_2>::computeDNDS(&xi[0],1, &DN(0,0));
      
      // compute jacobian
      matrix_type J = transpose(XX_)*DN;
      
      // compute function gradient
      vector_type gradF = transpose(J) * diff;
      for (UInt i=0; i<gradF.size(); ++i)
        grad[i] = gradF[i];
      
    }
    // return function value
    return 0.5 * transpose(diff)*diff;
  }
  
  point_type point() {
    vector_type x = (*this)(xi_);
    point_type p;
    for (UInt i=0; i<x.size(); ++i)
      p[i] = x[i];
    return p;
  }
  
  UInt iterations() const
  { return counter_; }

  
  static double wrap(const std::vector<double> &x, std::vector<double> &grad, void *data) {
    return (*reinterpret_cast<Distance_minimzer<_segment_2>*>(data))(x, grad); }
  
};

static Real constrain_triangle_3(const std::vector<Real> &xi, std::vector<Real> &grad, void *data) {
  if (!grad.empty()) {
    grad[0] = 1.;
    grad[1] = 1.;
  }
  return (xi[0] + xi[1] - 1.);
}

static void mf(unsigned m, double *result, unsigned n, const double *xi, double *grad, void *data) {
  
  assert(data == NULL);
  assert(n == 2);
  
  result[0] = -xi[0];
  result[1] = -xi[1];
  result[2] = xi[0] + xi[1] - 1.;
  if (grad) {
    grad[0] = grad[3] = -1.;
    grad[1] = grad[2] =  0.;
    grad[4] = grad[5] =  1.;
  }
}



template <>
class Distance_minimzer<_triangle_3> {
  
  static const UInt d = 3;
  static const UInt nb_nodes = 3;
  
  typedef Point<d> point_type;
  
  nlopt::opt opt_;         //!< Optimizator reference
  std::vector<Real> xi_;   //!< Master coordinate closest to point
  vector_type p_;          //!< Point to which the distance is minimized
  matrix_type XX_;         //!< Triangle coordinates
  UInt counter_;           //!< Optimization iteration counter
  Real fmin_;              //!< Minimum distance value
  
  void constructor_common() {
    
    // set optimization parameters
    std::vector<Real> lb(2); // default initialized to zero
    opt_.set_lower_bounds(lb);
    opt_.add_inequality_constraint(constrain_triangle_3, NULL, 1e-8);

    
//    std::vector<Real> tol(3,1e-8);
//    opt_.add_inequality_mconstraint(mf, NULL, tol);
    
    opt_.set_min_objective(wrap, this);
    opt_.set_ftol_abs(1e-8);
  }
  
public:
  
  //! Parameter constructor
  /*! \param r - Point coordinates
   * pts - Container of triangle points
   */
  template <class point_type, class point_container>
  Distance_minimzer(const point_type& p, const point_container& pts)
  : opt_(nlopt::LD_SLSQP, d-1), xi_(d-1), p_(d), XX_(nb_nodes,d), counter_() {
    
    // get triangle and point coordinates
    for (UInt i=0; i<d; ++i) {
      p_[i] = p[i];
      for (UInt j=0; j<nb_nodes; ++j)
        XX_(j,i) = pts[j][i];
    }

    // common constructor operation
    constructor_common();    
  }
  
  //! Parameter constructor
  /*! \param opt - Optimizer
   * \param r - Point coordinates
   * \param el - Finite element to which the distance is minimized
   * \param model - Solid mechanics model
   */
  Distance_minimzer(const Real *r, const Element *el, SolidMechanicsModel &model)
  : opt_(nlopt::LD_SLSQP, d-1), xi_(d-1), p_(d), XX_(nb_nodes,d), counter_() {
    
    // common constructor operation
    constructor_common();
    
    // get triangle coordinates from element and point coordinates
    Mesh& mesh = model.getMesh();
    const Array<Real> &X = model.getCurrentPosition();
    const Array<UInt> &conn = mesh.getConnectivity(el->type);
    for (UInt i=0; i<nb_nodes; ++i) {
      for (UInt j=0; j<d; ++j)
        XX_(i,j) = X(conn(el->element,i),j);
      p_(i) = r[i];
    }
  }
  
  void optimize() {
    
    // compute start point
    start();
        
#ifdef DEBUG_OPTIMIZE
    // optimize
    nlopt::result result = opt_.optimize(xi_, fmin_);
    if (result > 0)
      cout<<"Optimium found in "<<counter_<<" iterations: "<<fmin_<<endl;
    cout<<"Point at master coordinate "<<xi_[0]<<","<<xi_[1]<<": "<<point()<<endl;
    cout<<result<<endl;
#else
    // optimize
    opt_.optimize(xi_, fmin_);
#endif
  }
  
  point_type point() {
    vector_type x = (*this)(xi_);
    point_type p;
    for (UInt i=0; i<x.size(); ++i)
      p[i] = x[i];
    return p;
  }
  
  UInt iterations() const
  { return counter_; }
  
  const std::vector<Real>& master_coordinte()
  { return xi_; }

private:
  
  vector_type operator()(const std::vector<Real> &xi) {
    vector_type N(nb_nodes);
    ElementClass<_triangle_3>::computeShapes(&xi[0], &N(0));
    return transpose(XX_)*N;
  }
  
  void start() {
        
    Real min = std::numeric_limits<Real>::infinity();
    Real xstart[4][2] = { {0.,0.}, {1.,0.}, {0.,1.}, {1./3.,1./3.} }; // check center and corners of element
    int idx = -1;

    for (int i=0; i<4; ++i) {

      xi_[0] = xstart[i][0];
      xi_[1] = xstart[i][1];
      
      vector_type diff = (*this)(xi_) - p_;
      Real norm = diff.norm();
      if (norm < min) {
        min = norm;
        idx = i;
      }
    }
    xi_[0] = xstart[idx][idx];
    xi_[1] = xstart[idx][idx];
  }
  
  
  Real operator()(const std::vector<Real> &xi, std::vector<Real> &grad)
  {
    // increment function evaluation counter
    ++counter_;
    
    vector_type x = (*this)(xi);
    vector_type diff = x-p_;
    
    if (!grad.empty()) {
      
      // compute shape function derivatives
      matrix_type DN(nb_nodes, d-1);
      ElementClass<_triangle_3>::computeDNDS(&xi[0],&DN(0,0));
      
      // compute jacobian
      matrix_type J = transpose(XX_)*DN;
      
      // compute function gradient
      vector_type gradF = transpose(J) * diff;
      for (UInt i=0; i<gradF.size(); ++i)
        grad[i] = gradF[i];
    }
    // return function value
    return 0.5 * transpose(diff)*diff;
  }
  
  static double wrap(const std::vector<double> &x, std::vector<double> &grad, void *data) {
    return (*reinterpret_cast<Distance_minimzer<_triangle_3>*>(data))(x, grad); }
};




template <>
class Distance_minimzer<_triangle_6> {
  
  static const UInt d = 3;
  static const UInt nb_nodes = 6;
  
  typedef Point<d> point_type;
  
  nlopt::opt opt_;         //!< Optimizator reference
  std::vector<Real> xi_;   //!< Master coordinate closest to point
  vector_type p_;          //!< Point to which the distance is minimized
  matrix_type XX_;         //!< Triangle coordinates
  UInt counter_;           //!< Optimization iteration counter
  Real fmin_;              //!< Minimum distance value
  
  void constructor_common() {
    
//    std::vector<Real> tol(3,1e-2);
//    opt_.add_inequality_mconstraint(mf, NULL, tol);

    
    // set optimization parameters
    std::vector<Real> lb(2); // default initialized to zero
    opt_.set_lower_bounds(lb);
    // the quadratic triangle uses the same constraint as the linear triangle
    opt_.add_inequality_constraint(constrain_triangle_3, NULL, 1e-4);
    
    opt_.set_min_objective(wrap, this);

    opt_.set_ftol_abs(1e-2);
    opt_.set_ftol_rel(1e-2);
    opt_.set_stopval(Real());
    
    
    //    cout<<std::setprecision(10);
//    cout<<"\tstop val -> "<<opt_.get_stopval()<<endl;
//    cout<<"\tftol_rel -> "<<opt_.get_ftol_rel()<<endl;
//    cout<<"\tftol_abs -> "<<opt_.get_ftol_abs()<<endl;
//    cout<<"\txtol_rel -> "<<opt_.get_xtol_rel()<<endl;
//    cout<<"\tmaxeval -> "<<opt_.get_maxeval()<<endl;
//    cout<<"\tforce_stop -> "<<opt_.get_force_stop()<<endl;
//
//    
//    int w = 16;
//    
//    cout<<"#"<<std::setw(w-1)<<"x"<<std::setw(w)<<"y"<<std::setw(w)<<"xi"<<std::setw(w)<<"eta"<<std::setw(w)<<"x*"<<std::setw(w)<<"y*"<<std::setw(w)<<"f(xi)"<<std::setw(w)<<"iter"<<endl;
//    cout<<std::setprecision(3)<<std::fixed;

  }
  
public:
  
  //! Parameter constructor
  /*! \param r - Point coordinates
   * pts - Container of triangle points
   */
  template <class point_type, class point_container>
  Distance_minimzer(const point_type& p, const point_container& pts)
  : opt_(nlopt::LD_SLSQP, d-1), xi_(d-1), p_(d), XX_(nb_nodes,d), counter_() {

    // get triangle and point coordinates
    for (UInt i=0; i<d; ++i) {
      p_[i] = p[i];
      for (UInt j=0; j<nb_nodes; ++j)
        XX_(j,i) = pts[j][i];
    }

    // common constructor operation
    constructor_common();
  }
  
  //! Parameter constructor
  /*! \param opt - Optimizer
   * \param r - Point coordinates
   * \param el - Finite element to which the distance is minimized
   * \param model - Solid mechanics model
   */
  Distance_minimzer(const Real *r, const Element *el, SolidMechanicsModel &model)
  : opt_(nlopt::LD_SLSQP, d-1), xi_(d-1), p_(d), XX_(nb_nodes,d), counter_() {
    
    // common constructor operation
    constructor_common();
    
    // get triangle coordinates from element and point coordinates
    Mesh& mesh = model.getMesh();
    const Array<Real> &X = model.getCurrentPosition();
    const Array<UInt> &conn = mesh.getConnectivity(el->type);
    for (UInt i=0; i<nb_nodes; ++i) {
      for (UInt j=0; j<d; ++j)
        XX_(i,j) = X(conn(el->element,i),j);
      p_(i) = r[i];
    }
  }
  
  void optimize() {
    
    // compute start point
    start();
    
#ifdef DEBUG_OPTIMIZE
    // optimize
    nlopt::result result = opt_.optimize(xi_, fmin_);
    if (result > 0)
      cout<<"Optimium found in "<<counter_<<" iterations: "<<fmin_<<endl;
    cout<<"Point at master coordinate "<<xi_[0]<<","<<xi_[1]<<": "<<point()<<endl;
    cout<<result<<endl;
    if (result != 3)
      exit(1);
#else
    // optimize
    opt_.optimize(xi_, fmin_);
#endif
  }
  
  point_type point() {
    vector_type x = (*this)(xi_);
    point_type p;
    for (UInt i=0; i<x.size(); ++i)
      p[i] = x[i];
    return p;
  }
  
  UInt iterations() const
  { return counter_; }
  
  const std::vector<Real>& master_coordinte()
  { return xi_; }

private:
  
  vector_type operator()(const std::vector<Real> &xi) {
    vector_type N(nb_nodes);
    ElementClass<_triangle_6>::computeShapes(&xi[0], 1, &N(0));
    return transpose(XX_)*N;
  }
  
  void start() {
    
    Real min = std::numeric_limits<Real>::infinity();
    Real xstart[4][2] = { {0.,0.}, {1.,0.}, {0.,1.}, {1./3.,1./3.} }; // check center and corners of element
    int idx = -1;
    
    for (int i=0; i<4; ++i) {
      
      xi_[0] = xstart[i][0];
      xi_[1] = xstart[i][1];
      
      vector_type diff = (*this)(xi_) - p_;
      Real norm = diff.norm();
      if (norm < min) {
        min = norm;
        idx = i;
      }
    }
    xi_[0] = xstart[idx][idx];
    xi_[1] = xstart[idx][idx];
  }

    Real operator()(const std::vector<Real> &xi, std::vector<Real> &grad)
  {
    // increment function evaluation counter
    ++counter_;
        
    vector_type x = (*this)(xi);
    vector_type diff = x-p_;
    
    if (!grad.empty()) {
      
      // compute shape function derivatives
      matrix_type DN(nb_nodes, d-1);
      ElementClass<_triangle_6>::computeDNDS(&xi[0],1, &DN(0,0));
      
      // compute jacobian
      matrix_type J = transpose(XX_)*DN;
      
      // compute function gradient
      vector_type gradF = transpose(J) * diff;
      for (UInt i=0; i<gradF.size(); ++i)
        grad[i] = gradF[i];
    }
    
//    int w = 16;
//    cout<<std::setprecision(12);
//    cout<<std::setw(w)<<p_[0]<<std::setw(w)<<p_[1]<<std::setw(w)<<xi_[0]<<std::setw(w)<<xi_[1]<<std::setw(w)<<x[0]<<std::setw(w)<<x[1]<<std::setw(w)<<(0.5 * transpose(diff)*diff)<<std::setw(w)<<counter_<<endl;
    
    // return function value
    return 0.5 * transpose(diff)*diff;
  }
  
  static double wrap(const std::vector<double> &x, std::vector<double> &grad, void *data) {
    return (*reinterpret_cast<Distance_minimzer<_triangle_6>*>(data))(x, grad); }
};



__END_AKANTU__


#endif /* __AKANTU_OPTIMIZE_HH__ */


