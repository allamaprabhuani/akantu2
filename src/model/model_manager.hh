/**
 * @file   model_manager.hh
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @date   Wed Aug 22 12:14:00 2012
 *
 * @brief  higher order object that deals with collections of models
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

#ifndef __AKANTU_MODEL_MANAGER_HH__
#define __AKANTU_MODEL_MANAGER_HH__

#include <cstring>
#include <queue>



#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io_abaqus.hh"
#include "model.hh"
#include "aka_tree.hh"
#include "solid_mechanics_model.hh"
#include "aka_plane.hh"

#include "aka_timer.hh"

#define DEBUG_MANAGER 1


__BEGIN_AKANTU__


template <class Model_policy>
class Model_manager {
  
public:
  
  typedef Model_policy model_type;
  typedef model_type* model_pointer;
  typedef model_type& model_reference;
  
  typedef std::list<model_type*> model_container;
  typedef typename model_container::iterator model_iterator;
  typedef typename model_container::const_iterator const_model_iterator;
  
protected:
  
  model_container models_;            //!< Models
  
  
public:
  
  //! Default constructor
  Model_manager() : models_() {}
  
  virtual void add_model(model_reference m)
  { models_.push_back(&m); }
  
  virtual void add_model(model_pointer m)
  { models_.push_back(m); }
  
  friend std::ostream& operator<<(std::ostream& os, const Model_manager& mm) {
    
    os<<"Model manager info:"<<endl;
    os<<"  models: "<<mm.models_.size()<<endl;
    size_t i=0;
    for (const_model_iterator it = mm.models_.begin(); it != mm.models_.end(); ++it) {
      os<<"\tmodel "<<++i<<" memory address: "<<*it<<endl;
      os<<"\tmodel "<<**it<<endl;
    }
    
    return os;
  }
  
};

template <class tree>
void print_mathematica(const tree& t) {
  
  typedef typename tree::value_type volume_type;
  typedef typename tree::const_leaf_iterator iterator;
  
  int d = volume_type::dim();
  
  cout<<std::fixed;
  cout.precision(2);
  
  if (d == 2)
    cout<<"Graphics[{";
  else
    cout<<"Graphics3D[{";
  
  for (typename tree::const_leaf_iterator it = t.leaves_begin(); it != t.leaves_end(); ++it)
    cout<<*it<<", ";
  
  cout<<"}];"<<endl;
}


enum Kinematic_type { static_object_t = 0, dynamic_object_t = 1 };









template <class Bounding_policy, bool Accelerate = false, template <class> class Cost_policy = Cost_functor>
class Contact_model_manager : public Model_manager<SolidMechanicsModel> {
  
  typedef Real time_type;
  typedef ctimer chronograph_type;
  
  typedef SolidMechanicsModel model_type;
  typedef model_type* model_pointer;
  typedef model_type& model_reference;
  
  // Geometric
  typedef Bounding_policy volume_type;
  typedef typename volume_type::point_type point_type;
  typedef typename point_type::value_type value_type;
  typedef typename volume_type::aabb_type aabb_type;
  
  // Bounding volume hierarchy related types
  typedef Cost_policy<volume_type> cost_functor;
  typedef Tree<volume_type, Cost_policy > tree_type;
  typedef typename tree_type::leaf_iterator tree_leaf_iterator;
  typedef typename tree_type::iterator tree_iterator;
  typedef typename tree_type::const_iterator const_tree_iterator;
  typedef std::list<tree_type*> forest_container;
  typedef typename forest_container::iterator forest_iterator;
  typedef typename forest_container::const_iterator const_forest_iterator;
  
  // leaf information
  typedef std::vector<const Real*> coordinate_type;
  typedef std::map<tree_iterator, coordinate_type> leaves_data;
  typedef typename leaves_data::iterator leaves_data_iterator;
  typedef std::map<model_pointer, leaves_data> model_leaf_map;
  typedef typename model_leaf_map::iterator model_leaf_iterator;
  
  typedef std::tuple<leaves_data_iterator, leaves_data_iterator> check_type;
  
  struct Check_comparator {
    bool operator()(const check_type& c1, const check_type& c2) const
    { return std::get<0>(c1)->first < std::get<0>(c2)->first; }
  };
  
  typedef typename std::set<check_type, Check_comparator> check_container;
  typedef typename check_container::iterator check_iterator;
  
  
  // kinematic types
  typedef point_type velocity_type;
  typedef std::list<velocity_type> velocity_container;
  typedef typename velocity_container::iterator velocity_iterator;
  
  typedef unsigned long mask_size;
  
  // timer type
  typedef std::priority_queue<time_type, std::vector<time_type>, std::greater<time_type> > timer_type;
  
  //! Structure used to do a postorder update of tree hierarchies
  struct Updater {
    
    tree_type &t_;        //!< Reference to hierarchy
    
    Updater(tree_type& t) : t_(t) {}
    
    void operator()(tree_iterator it) {
      if (!it.is_leaf()) {
        volume_type& v = *it;
        volume_type& lv = *it.left();
        volume_type& rv = *it.right();
        v = lv + rv;
        assert(lv.last_time_ == rv.last_time_);
        v.last_time_ = lv.last_time_;
        v.velocity_ = 0.5 * (lv.velocity_ + rv.velocity_);
        v.acceleration_ = 0.5 * (lv.acceleration_ + rv.acceleration_);
      }
    }
  };
  
  struct Printer {
    void operator()(tree_iterator it)
    { cout<<*it<<", "; }
  };
  
  
  typedef std::tuple<time_type, tree_iterator, tree_iterator> tuple_type;
  
  struct Tuple_compare {
    bool operator()(const tuple_type& t1, const tuple_type& t2) const
    { return std::get<0>(t1) > std::get<0>(t2); }
  };
  
  typedef typename std::priority_queue<tuple_type, std::vector<tuple_type>, Tuple_compare> hierarchy_timer;
  
  
  template <bool flag>
  struct Bool2Type {
    enum { value = flag };
  };  
  
private:
  
  forest_container forest_;              //!< Bounding volume hierarchies
  hierarchy_timer timer_;                //!< Priority queue of times
  model_leaf_map leaf_data_;             //!< Leaf data
  mask_size masks_;                      //!< Variable used for static objects
  tree_iterator null_;
  check_container fine_check_;           //!< Container for fine check
  chronograph_type chrono_;
  
public:
  
  
  //! Default constructor
  Contact_model_manager()
  : Model_manager(), forest_(), timer_(), masks_(), null_(tree_iterator(nullptr)), chrono_() {
    timer_.push(std::make_tuple(time_type(), null_, null_));
    timer_.push(std::make_tuple(inf, null_, null_));
  }
  
  //! Destructor
  ~Contact_model_manager() {
    for (forest_iterator it = forest_.begin(); it != forest_.end(); ++it)
      delete *it;
  }
  
  virtual void add_model(model_pointer m, Kinematic_type k = dynamic_object_t) {
    
    m->initializeUpdateResidualData();
    models_.push_back(m);
    
    if (models_.size() > 8*sizeof(mask_size)) {
      cout<<"*** ERROR *** Type used for masks is too small to handle all models."<<endl;
      cout<<"Aborting..."<<endl;
      exit(1);
    }
    tree_type *tree = construct_tree_bottom_up<tree_type, volume_type>(*m, leaf_data_[m]);
    forest_.push_back(tree);
    
    //#ifdef DEBUG_MANAGER
    //    print_mathematica(*tree);
    //#endif
    
    // mask model as dynamic or static
    masks_ |= (k << (models_.size()-1));
  }
  
  
  void update_forest(time_type t) {
    
    cout<<"."<<endl;
    
    int k=0;
    model_leaf_iterator it = leaf_data_.begin();
    forest_iterator fit = forest_.begin();
    
    for (; it != leaf_data_.end(); ++it, ++fit) {
      
      // check if the object is dynamic to update
      if (!(masks_ & (1 << k++)))
        continue;
      
      // loop over leaves
      for (leaves_data_iterator lit = it->second.begin(); lit != it->second.end(); ++lit) {
        
        std::vector<point_type> pts;
        size_t nnodes = lit->second.size();
        pts.reserve(nnodes);
        for (size_t i=0; i<nnodes; ++i)
          pts.push_back(point_type(lit->second[i]));
        
        Real t_old = lit->first->last_time_;
        
        // get old position
        const point_type& p0 = lit->first->center();
        volume_type v = Volume_creator<volume_type>::create(lit->second);
        const point_type& p1 = v.center();
        
        const point_type& v0 = lit->first->velocity_;
        const point_type& v1 = v.velocity_;
        
        v.velocity_ = 1/(t - t_old) * (p1-p0);
        v.acceleration_ = 1/(t - t_old) * (v1-v0);
        
        v.last_time_ = t;
        
        // set new volume
        *lit->first = v;
      }
      
      tree_type &t = **fit;
      Updater u(t);
      postorder(t.root(),u);
    }
  }
  
  
  class Time_exception : public std::exception {
    virtual const char* what() const throw()
    { return "*** EXCEPTION *** Zero time increment."; }
  };
  
  
  /*! \param t - Current elapsed time
   * \param Dt - Time step
   */
  void intersect(time_type t, time_type& Dt) {
    
#ifdef DEBUG_MANAGER
    cout<<"-----------------------------------"<<endl;
    cout<<"Time ***"<<t<<endl;
#endif
    cout<<"Dt -> "<<Dt<<endl;
    
    static time_type Dt1 = Dt;
    
    const tuple_type& tuple = timer_.top();
    time_type top = std::get<0>(tuple);
    
    
    cout<<" Top time -> "<<top<<endl;
    
    Real offset = (Accelerate ? 2 : 1)*Dt;
    
    //    // update hierarchy
    //    if (t + offset > top && t > offset) {
    //
    //      // update hierarchies and get positions
    //      cout<<"updating BECAUSE t + offset > top"<<endl;
    //      update_forest(t);
    //
    //      cout<<"Graphics3D[{";
    //
    //      for (forest_iterator itt = forest_.begin(); itt != forest_.end(); ++itt)
    //        cout<<*(*itt)->root()<<",";
    //      cout<<"}]";
    //    }
    
    
    // update hierarchies and get positions
    update_forest(t);
    
    auto left = std::get<1>(tuple);
    auto right = std::get<2>(tuple);

    if (left != null_ && right != null_) {
    
    point_type cl = left->center();
    point_type cr = right->center();

      cout<<"sp1.Center=["<<cl[0]<<","<<cl[1]<<","<<cl[2]<<"]; sp1.Radius="<<left->radius()<<"; sp2.Center=["<<cr[0]<<","<<cr[1]<<","<<cr[2]<<"]; sp2.Radius="<<right->radius()<<";"<<endl;

      
    
//    cout<<"Sphere(); SetProperties(PhiResolution=32); SetProperties(ThetaResolution=32); SetProperties(Center=["<<cl[0]<<","<<cl[1]<<","<<cl[2]<<"]); SetProperties(Radius="<<left->radius()<<"); Show()"<<endl;
//
//    cout<<"Sphere(); SetProperties(PhiResolution=32); SetProperties(ThetaResolution=32); SetProperties(Center=["<<cr[0]<<","<<cr[1]<<","<<cr[2]<<"]); SetProperties(Radius="<<right->radius()<<"); Show()"<<endl;

    }
    
//    cout<<"ROOTS AFTER UPDATING"<<endl;
//    cout<<"Graphics3D[{Opacity[0.5],";
//    for (forest_iterator itt = forest_.begin(); itt != forest_.end(); ++itt)
//      cout<<*(*itt)->root()<<",";
//    cout<<"}]"<<endl;
    
    // set time step if required
    if (t + Dt > top && top != 0 && top != t) {
      
      Real tmp = Dt;
      Dt = std::abs(top - t);
      if (Dt > Dt1)
        Dt = Dt1;
      
      if (Dt < 1e-10) {
        cout<<"TOO SMALL"<<endl;
        throw Time_exception();
      }
      
      
      cout<<"*** SETTING Dt -> "<<Dt<<endl;
      
      //      if (Dt == 0) {
      //
      //        // update hierarchies and get positions
      //        update_forest(t);
      //
      //        cout<<"Graphics3D[{";
      //        for (forest_iterator itt = forest_.begin(); itt != forest_.end(); ++itt)
      //          cout<<*(*itt)->root()<<",";
      //        cout<<"}]"<<endl;
      //        throw Time_exception();
      //
      //      }
      
      cout<<"new Dt -> "<<Dt<<endl;
      for (model_iterator mit = models_.begin(); mit != models_.end(); ++mit) {
        (*mit)->setTimeStep(Dt);
        cout<<"calling to set time step "<<(Dt)<<endl;
      }
    }
    
    // if there are enough models to carry out intersection
    if (models_.size() > 1 && top <= t) {
      
      //      // update hierarchies and get positions
      //      update_forest(t);
      
      
      // if first step, add next time to timer and return because there
      // is not enough information for the computation of intersection times
      if (t <= (Accelerate ? 2 : 1)*Dt) {
        // remove time from timer
        timer_.pop();
        timer_.push(std::make_tuple(t, null_, null_));
        cout<<chrono_<<endl;
        return;
      }
      
      
#ifdef DEBUG_MANAGER
      cout<<"top time -> "<<top<<endl;
      cout<<"timer size -> "<<timer_.size()<<endl;
#endif
      
      // get collision time
      //      t = std::min(t, top);
      
      // check if iterators are null to compute O(n^2) collision times
      tree_iterator it1 = std::get<1>(tuple);
      tree_iterator it2 = std::get<2>(tuple);
      
      // remove time from timer
      timer_.pop();
      
      if (it1 == null_ || it2 == null_) {
        
        // do O(n^2) operation to obtain next time of intersections
        for (forest_iterator it1 = forest_.begin(); it1 != --forest_.end(); ++it1) {
          
          forest_iterator it2 = it1;
          for (++it2; it2 != forest_.end(); ++it2) {
            
            // get collision time
#ifdef DEBUG_MANAGER
            cout<<"Calling check_collision in O(n^2) branch"<<endl;
#endif
            check_collision(t, Dt, (*it1)->root(), (*it2)->root());
            
          } // inner hierarchy loop
        } // outer hierarchy loop
        
      }
      // else use collision information previously computed (avoids O(n^2) operation above)
      else {
#ifdef DEBUG_MANAGER
        cout<<"Calling check_collision in non-O(n^2) branch"<<endl;
#endif
        check_collision(t, Dt, it1, it2);
      }
      
    } // if statement on enough models
    // else do nothing as there are not enough models to carry out intersection
    // or the check engine is shut down until the next time in timer
    
    // carry out fine check if necessary
    if (!fine_check_.empty()) {
      
      // loop over elements in the set
      for (check_iterator it = fine_check_.begin(); it != fine_check_.end(); ++it) {
        
        const check_type& c = *it;
        leaves_data_iterator lit1 = std::get<0>(c);
        leaves_data_iterator lit2 = std::get<1>(c);
        
        bool penetrate = check_penetration(lit1, lit2);
        if (penetrate) {
          cout<<"FINALLY PENETRATING at time " <<t<<" !!!"<<endl;
          exit(1);
        }
      }
    } // non-empty fine check container
    
    
    const tuple_type& tuple2 = timer_.top();
    time_type top2 = std::get<0>(tuple);
    
    cout<<"top -> "<<top<<endl;
    cout<<"top2 -> "<<top2<<endl;
    
    if (top2 - top > 0) {
      Dt = std::min(std::abs(top2 - top) + 1e-6, Dt1);
      cout<<"***** CHANGING DT TO "<<Dt<<endl;
      for (model_iterator mit = models_.begin(); mit != models_.end(); ++mit) {
        (*mit)->setTimeStep(Dt);
        cout<<"calling to set time step "<<(Dt)<<endl;
      }
    }
  }
  
  
  
  // checks if a point s moves away from a plane, indicating that a node has
  // penetrated and that further checking is required
  bool check_penetration(const point_type& s, const point_type& v, Plane& p) {
    
    // compute distance to plane
    value_type dist = p.normal() * s - p.distance();
    value_type denom = p.normal() * v;
    
    // point is moving away from the plane, indicating that has penetrated already
    if (denom * dist > 0.)
      return true;
    return false;
  }
  
  
private:
  
  template <class iterator>
  void check_collision(time_type t, time_type& Dt, iterator it1, iterator it2) {
    
    check_collision(t, Dt, it1, it2, Bool2Type<Accelerate>());
    
    while (!timer_.empty()) {
      
      const tuple_type& tuple = timer_.top();
      
      time_type tstar = std::get<0>(tuple);
      iterator left = std::get<1>(tuple);
      iterator right = std::get<2>(tuple);
      
      // return if collision happens later in time
      if (tstar > t)
        return;
      
      check_collision(t, Dt, left, right, Bool2Type<Accelerate>());
      
      // remove tuple from the timer and call chack collision on their
      // iterators
      timer_.pop();
    }
  }
  
  
  template <class iterator>
  void handle_collision(time_type t, time_type& Dt, iterator it1, iterator it2) {
    
    cout<<" INSIDE HANDLE COLLIISON"<<endl;
    
    // if volumes are leaves
    if (it1.is_leaf() && it2.is_leaf()) {
      
#ifdef DEBUG_MANAGER
      cout<<"    Leaves found"<<endl;
    
      const point_type& cl = it1->center();
      const point_type& cr = it2->center();
            
      cout<<"sp1.Center=["<<cl[0]<<","<<cl[1]<<","<<cl[2]<<"]; sp1.Radius="<<it1->radius()<<"; sp2.Center=["<<cr[0]<<","<<cr[1]<<","<<cr[2]<<"]; sp2.Radius="<<it2->radius()<<";"<<endl;

#endif
      
      
      

      // add leaves to carry out penetration tests
      leaves_data_iterator lit1, lit2;
      
      // check for penetration
      for (model_leaf_iterator it = leaf_data_.begin(); it != leaf_data_.end(); ++it) {
        
        leaves_data_iterator lit = it->second.find(it1);
        if (lit != it->second.end())
          lit1 = lit;
        
        lit = it->second.find(it2);
        if (lit != it->second.end())
          lit2 = lit;
      }
      
      fine_check_.insert(std::make_tuple(lit1,lit2));
      return;
      
      //        // UNCOMMENT TO WORK WITH AXIS-ALIGNED BOUNDING BOXES
      //        // now that leaves have been found, check collision of axis-aligned bounding boxes
      //        aabb_type bb1, bb2;
      //
      //        for (model_leaf_iterator it = leaf_data_.begin(); it != leaf_data_.end(); ++it) {
      //
      //          leaves_data_iterator lit = it->second.find(it1);
      //          if (lit != it->second.end())
      //            for (int i=0; i<lit->second.size(); ++i)
      //              bb1 += point_type(lit->second[i]);
      //
      //          lit = it->second.find(it2);
      //          if (lit != it->second.end())
      //            for (int i=0; i<lit->second.size(); ++i)
      //              bb2 += point_type(lit->second[i]);
      //        }
      //
      //        // assign velocities
      //        bb1.velocity_ = s1.velocity_;
      //        bb2.velocity_ = s2.velocity_;
      //
      //        // get time of intersection between aabbs
      //        Real ts = check_collision(bb1,bb2);
    }
    
    // traverse down the left hierarchy if its volume is bigger or if the
    // right volume is a leaf
    if ((!it1.is_leaf() && it1->measure() >= it2->measure()) || it2.is_leaf()) {
      
#ifdef DEBUG_MANAGER
      cout<<"    s1 >= s2 OR s2 is leaf"<<endl;
#endif
      iterator lit = it1.left();
      iterator rit = it1.right();
      assert(lit != null_);
      assert(rit != null_);
      check_collision(t, Dt, lit, it2, Bool2Type<Accelerate>());
      check_collision(t, Dt, rit, it2, Bool2Type<Accelerate>());
    }
    // traverse down the right hierarchy if its volume is bigger or if the
    // right volume is a leaf
    else if ((!it2.is_leaf() && it2->measure() > it1->measure()) || it1.is_leaf()) {
      
#ifdef DEBUG_MANAGER
      cout<<"    s1 < s2 OR s1 is leaf"<<endl;
#endif
      
      iterator lit = it2.left();
      iterator rit = it2.right();
      assert(lit != null_);
      assert(rit != null_);
      check_collision(t, Dt, it1, lit, Bool2Type<Accelerate>());
      check_collision(t, Dt, it1, rit, Bool2Type<Accelerate>());
    }
    else {
      // should never call this branch
      assert (false);
    }
  }
  
  
  
  // check collision considering acceleration
  // the equation to solve is
  //
  //   alpha t^4/4 + beta t^3 + gamma t^2 + delta t + epsilon = 0
  //
  // where alpha = (a.a)/4, beta = (a.v), gamma = a.s + v.v, delta = 2(s.v),
  // epsilon = s.s - r^2, a = a1-a2, v = v1-v2, s = c1-c2, r = r1-r2
  //
  template <class iterator>
  void check_collision(time_type t, time_type& Dt, iterator it1, iterator it2, Bool2Type<true>) {
    
    
    // check for collision
    if (*it1 & *it2)
      handle_collision(t, Dt, it1, it2);
    
    const point_type& v1 = it1->velocity_;
    const point_type& v2 = it2->velocity_;
    
    point_type a1 = it1->acceleration_;
    point_type a2 = it2->acceleration_;
    
    // TODO REMOVE HARDCODED ACCELERATION
    a1[0] = a1[1] = a2[0] = a2[1] = 0;
    a1[2] = a1[2] = -9.81;
    
    Real r = it1->radius() + it2->radius();
    point_type s = it2->center() - it1->center();
    point_type v = v2 - v1;
    point_type a = a2 - a1;
    
    // get coefficients
    Real alpha = (a*a)/4;
    Real beta = a*v;
    Real gamma = a*s + v*v;
    Real delta = 2*(s*v);
    Real epsilon = s*s - r*r;
    
    // obtain roots from quartic equation by calling the function such that
    // the coefficient for the quartic term is equal to one
    std::vector<Real> x(4, -1);
    
    cout<<"calling solve quartic with coefficients "<<endl;
    cout<<(beta/alpha)<<", "<<(gamma/alpha)<<", "<<(delta/alpha)<<", "<<(epsilon/alpha)<<endl;
    uint32_t roots = solve_quartic(beta/alpha, gamma/alpha, delta/alpha, epsilon/alpha,
                                   &x[0], &x[1], &x[2], &x[3]);
    
    // if there are roots, take the first one as an indication of the collision
    if (roots > 0) {
      
      Real tmin = std::numeric_limits<Real>::infinity();
      for (size_t i=0; i<x.size(); ++i)
        if (x[i] > 0)
          tmin = std::min(tmin, x[i]);
      if (tmin+t != std::numeric_limits<Real>::quiet_NaN()) {
        timer_.push(std::make_tuple(tmin+t, it1, it2));
        // change time step
        cout<<"APPROXIMATE COLLISION TIME -> t + tmin = "<<t<<" + "<<tmin<<" = "<<(tmin+t)<<endl;
        if (tmin == 0) {
          assert(false);
        }
        
        //        if (tmin < Dt) {
        //          cout<<"CHANGING TIME STEP TO "<<tmin<<endl;
        //          Dt = tmin;
        //        }
      }
    }
    
    cout<<" roots -> "<<roots<<endl;
    cout<<" values -> ";
    for (size_t i=0; i<x.size(); ++i)
      cout<<x[i]<<", ";
    cout<<endl;
  }
  
  // check collision for constant velocity
  template <class iterator>
  void check_collision(time_type t, time_type& Dt, iterator it1, iterator it2, Bool2Type<false> flag) {
    
    const volume_type& s1 = *it1;
    const volume_type& s2 = *it2;
    
    const point_type& v1 = s1.velocity_;
    const point_type& v2 = s2.velocity_;
    
    // vector between spheress
    point_type s = s2.center() - s1.center();
    // relative motion of s1 with respect to stationary s0
    point_type v = v2 - v1;
    // sum of radii
    value_type r = s1.radius() + s2.radius();
    value_type c = s*s - r*r;
    
#ifdef DEBUG_MANAGER
    cout<<"Checking collisiong between:"<<endl;
    cout<<"  "<<s1<<", velocity: "<<v1<<endl;
    cout<<"  "<<s2<<", velocity: "<<v2<<endl;
    cout<<"  Relative velocity: "<<v<<endl;
#endif
    
    // already intersecting, call recursively on branches
    if (c < 0.) {
      
#ifdef DEBUG_MANAGER
      cout<<"  Intersecting bounding volumes"<<endl;
#endif
      handle_collision(t, Dt, it1, it2);
      
    }
    
    value_type epsilon = 1e-8;
    value_type a = v*v;
    // if spheres not moving relative to each other
    if (a < epsilon) {
#ifdef DEBUG_MANAGER
      cout<<"    Objects not moving relative to each other"<<endl;
      cout<<"    a: "<<a<<endl;
#endif
      return;
    }
    
    value_type b = v*s;
    // if spheres not moving towards each other
    if (b >= 0.) {
#ifdef DEBUG_MANAGER
      cout<<"    Objects not moving towards each other"<<endl;
      cout<<"    b: "<<b<<endl;
#endif
      return;
    }
    
    value_type d = b*b - a*c;
    // if no real-valued root (d < 0), spheres do not intersect
    // otherwise add time to timer
    if (d >= 0.) {
      time_type ts = (-b - sqrt(d))/a;
      if (ts > 0.) {
#ifdef DEBUG_MANAGER
        cout<<"    Objects intersects at time "<<(t+ts)<<endl;
#endif
        timer_.push(std::make_tuple(t + ts, it1, it2));
      }
#ifdef DEBUG_MANAGER
      else {
        cout<<"    ts negative: "<<ts<<endl;
      }
#endif
      
    }
#ifdef DEBUG_MANAGER
    else {
      cout<<"    Objects do not intersect"<<endl;
      cout<<"    d: "<<d<<endl;
    }
#endif
  }
  
  bool check_penetration(leaves_data_iterator it1, leaves_data_iterator it2) {
    
    const coordinate_type& c1 = it1->second;
    const coordinate_type& c2 = it2->second;
    
    for (size_t i=0; i<c1.size(); ++i) {
      
      point_type r(c1[i]);
      
      // swith on the number of coordinates of the second container
      switch (c2.size()) {
          
          // triangle case
        case 3:
        {
          // form plane from three points
          
          // get element extreme points
          point_type o(c2[0]);
          point_type p(c2[1]);
          point_type q(c2[2]);
          
          Plane pi(o,p,q);
          
          // check for actual penetration using relative velocity
          if (check_penetration(r, it1->first->velocity_ - it2->first->velocity_, pi))
            return true;
        }
          break;
          
        default:
          cout<<"*** ERROR *** Function check_penetration in file "<<__FILE__<<", line "<<__LINE__;
          cout<<", is not implemented for node container of size "<<c2.size()<<endl;
          cout<<"*** ABORTING *** "<<endl;
          exit(1);
          break;
      }
    } // loop over coordinates of first container
    
    return false;
  }
  
  // intersects aabbs a and b moving with constant velocities.
  // On intersection, return time of first contact
  time_type check_collision(const aabb_type& a, const aabb_type& b) {
    
    // overlapping aabbs
    if (a & b) {
      cout<<"*** ERROR *** Colliding bounding boxes"<<endl;
      cout<<"              Aborting..."<<endl;
      exit(1);
    }
    
    // compute relative velocity
    const point_type v = b.velocity_ - a.velocity_;
    
    // initialize times of first and last contact
    time_type tf = 0., tl = 0.;
    
    // for each axis, determine the time of first and last contact
    for (int i=0; i<point_type::dim(); ++i) {
      
      if (v[i] < 0.) {
        if (b.max(i) < a.min(i))
          return inf; // non-intersecting
        if (a.max(i) < b.min(i))
          tf = std::max((a.max(i) - b.min(i)) / v[i], tf);
        if (b.max(i) > a.min(i))
          tl = std::min((a.min(i) - b.max(i)) / v[i], tl);
      }
      if (v[i] > 0.) {
        if (b.min(i) > a.max(i))
          return inf; // non-intersecting
        if (b.max(i) < a.min(i))
          tf = std::max((a.min(i) - b.max(i)) / v[i], tf);
        if (a.max(i) > b.min(i))
          tl = std::min((a.max(i) - b.min(i)) / v[i], tl);
      }
    }
    return tf;
  }
  
  /*! \brief Solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d = 0.
   *
   *  Solves the quartic equation. Returns the number of roots
   *  found. Roots are filled in ascending order.
   *  Code taken from the Ion Beam Simulator, which is distrubuted under
   *  the terms of the GNU General Public License as published by the Free
   *  Software Foundation; either version 2 of the License, or (at your
   *  option) any later version.
   */
  uint32_t solve_quartic(Real a, Real b, Real c, Real d,
                         Real *x0, Real *x1, Real *x2, Real *x3 );
  
  
  friend std::ostream& operator<<(std::ostream& os, const Contact_model_manager& mm) {
    
    os<<"Contact model manager info:"<<endl;
    os<<"  models: "<<mm.models_.size()<<endl;
    size_t i=0;
    const_forest_iterator tit = mm.forest_.begin();
    for (const_model_iterator it = mm.models_.begin(); it != mm.models_.end(); ++it) {
      os<<"\tmodel "<<++i<<" memory address: "<<*it<<endl;
      os<<"\tmodel: "<<**it<<endl;
      os<<"\ttree: ";
      print_mathematica(**tit++);
    }
    return os;
  }
};





// NEW IMPLEMENTATION
template <class Bounding_policy, bool Accelerate = false, template <class> class Cost_policy = Cost_functor>
class Contact_model_manager2 : public Model_manager<SolidMechanicsModel> {
  
  typedef Real time_type;
  typedef ctimer chronograph_type;
  
  typedef SolidMechanicsModel model_type;
  typedef model_type* model_pointer;
  typedef model_type& model_reference;
  
  // Geometric
  typedef Bounding_policy volume_type;
  typedef typename volume_type::point_type point_type;
  typedef typename point_type::value_type value_type;
  typedef typename volume_type::aabb_type aabb_type;
  
  // Bounding volume hierarchy related types
  typedef Cost_policy<volume_type> cost_functor;
  typedef Tree<volume_type, Cost_policy > tree_type;
  typedef typename tree_type::leaf_iterator tree_leaf_iterator;
  typedef typename tree_type::iterator tree_iterator;
  typedef typename tree_type::const_iterator const_tree_iterator;
  typedef std::list<tree_type*> forest_container;
  typedef typename forest_container::iterator forest_iterator;
  typedef typename forest_container::const_iterator const_forest_iterator;
  
  // leaf information
  typedef std::vector<const Real*> coordinate_type;
  typedef std::map<tree_iterator, coordinate_type> leaves_data;
  typedef typename leaves_data::iterator leaves_data_iterator;
  typedef std::map<model_pointer, leaves_data> model_leaf_map;
  typedef typename model_leaf_map::iterator model_leaf_iterator;
  
  typedef std::tuple<leaves_data_iterator, leaves_data_iterator> check_type;
  
  struct Check_comparator {
    bool operator()(const check_type& c1, const check_type& c2) const
    { return std::get<0>(c1)->first < std::get<0>(c2)->first; }
  };
  
  typedef typename std::set<check_type, Check_comparator> check_container;
  typedef typename check_container::iterator check_iterator;
  
  
  // kinematic types
  typedef point_type velocity_type;
  typedef std::list<velocity_type> velocity_container;
  typedef typename velocity_container::iterator velocity_iterator;
  
  typedef unsigned long mask_size;
  
  // timer type
  typedef std::priority_queue<time_type, std::vector<time_type>, std::greater<time_type> > timer_type;
  
  //! Structure used to do a postorder update of tree hierarchies
  struct Updater {
    
    tree_type &t_;        //!< Reference to hierarchy
    
    Updater(tree_type& t) : t_(t) {}
    
    void operator()(tree_iterator it) {
      if (!it.is_leaf()) {
        volume_type& v = *it;
        volume_type& lv = *it.left();
        volume_type& rv = *it.right();
        v = lv + rv;
        assert(lv.last_time_ == rv.last_time_);
        v.last_time_ = lv.last_time_;
        v.velocity_ = 0.5 * (lv.velocity_ + rv.velocity_);
        v.acceleration_ = 0.5 * (lv.acceleration_ + rv.acceleration_);
      }
    }
  };
  
  struct Printer {
    void operator()(tree_iterator it)
    { cout<<*it<<", "; }
  };
  
  
  typedef std::tuple<time_type, tree_iterator, tree_iterator> tuple_type;
  
  struct Tuple_compare {
    bool operator()(const tuple_type& t1, const tuple_type& t2) const
    { return std::get<0>(t1) > std::get<0>(t2); }
  };
  
  typedef typename std::priority_queue<tuple_type, std::vector<tuple_type>, Tuple_compare> hierarchy_timer;
  
  
  template <bool flag>
  struct Bool2Type {
    enum { value = flag };
  };
  
private:
  
  forest_container forest_;              //!< Bounding volume hierarchies
  hierarchy_timer timer_;                //!< Priority queue of times
  model_leaf_map leaf_data_;             //!< Leaf data
  mask_size masks_;                      //!< Variable used for static objects
  tree_iterator null_;
  check_container fine_check_;           //!< Container for fine check
  chronograph_type chrono_;
  
public:
  
  
  //! Default constructor
  Contact_model_manager2()
  : Model_manager(), forest_(), timer_(), masks_(), null_(tree_iterator(nullptr)), chrono_() {
    timer_.push(std::make_tuple(time_type(), null_, null_));
    timer_.push(std::make_tuple(inf, null_, null_));
  }
  
  //! Destructor
  ~Contact_model_manager2() {
    for (forest_iterator it = forest_.begin(); it != forest_.end(); ++it)
      delete *it;
  }
  
  virtual void add_model(model_pointer m, Kinematic_type k = dynamic_object_t) {
    
    m->initializeUpdateResidualData();
    models_.push_back(m);
    
    if (models_.size() > 8*sizeof(mask_size)) {
      cout<<"*** ERROR *** Type used for masks is too small to handle all models."<<endl;
      cout<<"Aborting..."<<endl;
      exit(1);
    }
    tree_type *tree = construct_tree_bottom_up<tree_type, volume_type>(*m, leaf_data_[m]);
    forest_.push_back(tree);
    
    //#ifdef DEBUG_MANAGER
    //    print_mathematica(*tree);
    //#endif
    
    // mask model as dynamic or static
    masks_ |= (k << (models_.size()-1));
  }
  
  
  void update_forest(time_type t) {
    
    cout<<"."<<endl;
    
    int k=0;
    model_leaf_iterator it = leaf_data_.begin();
    forest_iterator fit = forest_.begin();
    
    for (; it != leaf_data_.end(); ++it, ++fit) {
      
      // check if the object is dynamic to update
      if (!(masks_ & (1 << k++)))
        continue;
      
      // loop over leaves
      for (leaves_data_iterator lit = it->second.begin(); lit != it->second.end(); ++lit) {
        
        std::vector<point_type> pts;
        size_t nnodes = lit->second.size();
        pts.reserve(nnodes);
        for (size_t i=0; i<nnodes; ++i)
          pts.push_back(point_type(lit->second[i]));
        
        Real t_old = lit->first->last_time_;
        
        // get old position
        const point_type& p0 = lit->first->center();
        volume_type v = Volume_creator<volume_type>::create(lit->second);
        const point_type& p1 = v.center();
        
        const point_type& v0 = lit->first->velocity_;
        const point_type& v1 = v.velocity_;
        
        v.velocity_ = 1/(t - t_old) * (p1-p0);
        v.acceleration_ = 1/(t - t_old) * (v1-v0);
        
        v.last_time_ = t;
        
        // set new volume
        *lit->first = v;
      }
      
      tree_type &t = **fit;
      Updater u(t);
      postorder(t.root(),u);
    }
  }
  
  
  class Time_exception : public std::exception {
    virtual const char* what() const throw()
    { return "*** EXCEPTION *** Zero time increment."; }
  };
  
  
  /*! \param t - Current elapsed time
   * \param Dt - Time step
   */
  void intersect(time_type t, time_type& Dt) {
    
#ifdef DEBUG_MANAGER
    cout<<"-----------------------------------"<<endl;
    cout<<"Time ***"<<t<<endl;
#endif
    cout<<"Dt -> "<<Dt<<endl;
    
    static time_type Dt1 = Dt;
    
    const tuple_type& tuple = timer_.top();
    time_type top = std::get<0>(tuple);
    
    
    cout<<" Top time -> "<<top<<endl;
    
    Real offset = (Accelerate ? 2 : 1)*Dt;
    
    //    // update hierarchy
    //    if (t + offset > top && t > offset) {
    //
    //      // update hierarchies and get positions
    //      cout<<"updating BECAUSE t + offset > top"<<endl;
    //      update_forest(t);
    //
    //      cout<<"Graphics3D[{";
    //
    //      for (forest_iterator itt = forest_.begin(); itt != forest_.end(); ++itt)
    //        cout<<*(*itt)->root()<<",";
    //      cout<<"}]";
    //    }
    
    
    // update hierarchies and get positions
    update_forest(t);
    
    auto left = std::get<1>(tuple);
    auto right = std::get<2>(tuple);
    
    if (left != null_ && right != null_) {
      
      point_type cl = left->center();
      point_type cr = right->center();
      
      cout<<"sp1.Center=["<<cl[0]<<","<<cl[1]<<","<<cl[2]<<"]; sp1.Radius="<<left->radius()<<"; sp2.Center=["<<cr[0]<<","<<cr[1]<<","<<cr[2]<<"]; sp2.Radius="<<right->radius()<<";"<<endl;
      
      
      
      //    cout<<"Sphere(); SetProperties(PhiResolution=32); SetProperties(ThetaResolution=32); SetProperties(Center=["<<cl[0]<<","<<cl[1]<<","<<cl[2]<<"]); SetProperties(Radius="<<left->radius()<<"); Show()"<<endl;
      //
      //    cout<<"Sphere(); SetProperties(PhiResolution=32); SetProperties(ThetaResolution=32); SetProperties(Center=["<<cr[0]<<","<<cr[1]<<","<<cr[2]<<"]); SetProperties(Radius="<<right->radius()<<"); Show()"<<endl;
      
    }
    
    //    cout<<"ROOTS AFTER UPDATING"<<endl;
    //    cout<<"Graphics3D[{Opacity[0.5],";
    //    for (forest_iterator itt = forest_.begin(); itt != forest_.end(); ++itt)
    //      cout<<*(*itt)->root()<<",";
    //    cout<<"}]"<<endl;
    
    // set time step if required
    if (t + Dt > top && top != 0 && top != t) {
      
      Real tmp = Dt;
      Dt = std::abs(top - t);
      if (Dt > Dt1)
        Dt = Dt1;
      
      if (Dt < 1e-10) {
        cout<<"TOO SMALL"<<endl;
        throw Time_exception();
      }
      
      
      cout<<"*** SETTING Dt -> "<<Dt<<endl;
      
      //      if (Dt == 0) {
      //
      //        // update hierarchies and get positions
      //        update_forest(t);
      //
      //        cout<<"Graphics3D[{";
      //        for (forest_iterator itt = forest_.begin(); itt != forest_.end(); ++itt)
      //          cout<<*(*itt)->root()<<",";
      //        cout<<"}]"<<endl;
      //        throw Time_exception();
      //
      //      }
      
      cout<<"new Dt -> "<<Dt<<endl;
      for (model_iterator mit = models_.begin(); mit != models_.end(); ++mit) {
        (*mit)->setTimeStep(Dt);
        cout<<"calling to set time step "<<(Dt)<<endl;
      }
    }
    
    // if there are enough models to carry out intersection
    if (models_.size() > 1 && top <= t) {
      
      //      // update hierarchies and get positions
      //      update_forest(t);
      
      
      // if first step, add next time to timer and return because there
      // is not enough information for the computation of intersection times
      if (t <= (Accelerate ? 2 : 1)*Dt) {
        // remove time from timer
        timer_.pop();
        timer_.push(std::make_tuple(t, null_, null_));
        cout<<chrono_<<endl;
        return;
      }
      
      
#ifdef DEBUG_MANAGER
      cout<<"top time -> "<<top<<endl;
      cout<<"timer size -> "<<timer_.size()<<endl;
#endif
      
      // get collision time
      //      t = std::min(t, top);
      
      // check if iterators are null to compute O(n^2) collision times
      tree_iterator it1 = std::get<1>(tuple);
      tree_iterator it2 = std::get<2>(tuple);
      
      // remove time from timer
      timer_.pop();
      
      if (it1 == null_ || it2 == null_) {
        
        // do O(n^2) operation to obtain next time of intersections
        for (forest_iterator it1 = forest_.begin(); it1 != --forest_.end(); ++it1) {
          
          forest_iterator it2 = it1;
          for (++it2; it2 != forest_.end(); ++it2) {
            
            // get collision time
#ifdef DEBUG_MANAGER
            cout<<"Calling check_collision in O(n^2) branch"<<endl;
#endif
            check_collision(t, Dt, (*it1)->root(), (*it2)->root());
            
          } // inner hierarchy loop
        } // outer hierarchy loop
        
      }
      // else use collision information previously computed (avoids O(n^2) operation above)
      else {
#ifdef DEBUG_MANAGER
        cout<<"Calling check_collision in non-O(n^2) branch"<<endl;
#endif
        check_collision(t, Dt, it1, it2);
      }
      
    } // if statement on enough models
    // else do nothing as there are not enough models to carry out intersection
    // or the check engine is shut down until the next time in timer
    
    // carry out fine check if necessary
    if (!fine_check_.empty()) {
      
      // loop over elements in the set
      for (check_iterator it = fine_check_.begin(); it != fine_check_.end(); ++it) {
        
        const check_type& c = *it;
        leaves_data_iterator lit1 = std::get<0>(c);
        leaves_data_iterator lit2 = std::get<1>(c);
        
        bool penetrate = check_penetration(lit1, lit2);
        if (penetrate) {
          cout<<"FINALLY PENETRATING at time " <<t<<" !!!"<<endl;
          exit(1);
        }
      }
    } // non-empty fine check container
    
    
    const tuple_type& tuple2 = timer_.top();
    time_type top2 = std::get<0>(tuple);
    
    cout<<"top -> "<<top<<endl;
    cout<<"top2 -> "<<top2<<endl;
    
    if (top2 - top > 0) {
      Dt = std::min(std::abs(top2 - top) + 1e-6, Dt1);
      cout<<"***** CHANGING DT TO "<<Dt<<endl;
      for (model_iterator mit = models_.begin(); mit != models_.end(); ++mit) {
        (*mit)->setTimeStep(Dt);
        cout<<"calling to set time step "<<(Dt)<<endl;
      }
    }
  }
  
  
  
  // checks if a point s moves away from a plane, indicating that a node has
  // penetrated and that further checking is required
  bool check_penetration(const point_type& s, const point_type& v, Plane& p) {
    
    // compute distance to plane
    value_type dist = p.normal() * s - p.distance();
    value_type denom = p.normal() * v;
    
    // point is moving away from the plane, indicating that has penetrated already
    if (denom * dist > 0.)
      return true;
    return false;
  }
  
  
private:
  
  template <class iterator>
  void check_collision(time_type t, time_type& Dt, iterator it1, iterator it2) {
    
    check_collision(t, Dt, it1, it2, Bool2Type<Accelerate>());
    
    while (!timer_.empty()) {
      
      const tuple_type& tuple = timer_.top();
      
      time_type tstar = std::get<0>(tuple);
      iterator left = std::get<1>(tuple);
      iterator right = std::get<2>(tuple);
      
      // return if collision happens later in time
      if (tstar > t)
        return;
      
      check_collision(t, Dt, left, right, Bool2Type<Accelerate>());
      
      // remove tuple from the timer and call chack collision on their
      // iterators
      timer_.pop();
    }
  }
  
  
  template <class iterator>
  void handle_collision(time_type t, time_type& Dt, iterator it1, iterator it2) {
    
    cout<<" INSIDE HANDLE COLLIISON"<<endl;
    
    // if volumes are leaves
    if (it1.is_leaf() && it2.is_leaf()) {
      
#ifdef DEBUG_MANAGER
      cout<<"    Leaves found"<<endl;
      
      const point_type& cl = it1->center();
      const point_type& cr = it2->center();
      
      cout<<"sp1.Center=["<<cl[0]<<","<<cl[1]<<","<<cl[2]<<"]; sp1.Radius="<<it1->radius()<<"; sp2.Center=["<<cr[0]<<","<<cr[1]<<","<<cr[2]<<"]; sp2.Radius="<<it2->radius()<<";"<<endl;
      
#endif
      
      
      
      
      // add leaves to carry out penetration tests
      leaves_data_iterator lit1, lit2;
      
      // check for penetration
      for (model_leaf_iterator it = leaf_data_.begin(); it != leaf_data_.end(); ++it) {
        
        leaves_data_iterator lit = it->second.find(it1);
        if (lit != it->second.end())
          lit1 = lit;
        
        lit = it->second.find(it2);
        if (lit != it->second.end())
          lit2 = lit;
      }
      
      fine_check_.insert(std::make_tuple(lit1,lit2));
      return;
      
      //        // UNCOMMENT TO WORK WITH AXIS-ALIGNED BOUNDING BOXES
      //        // now that leaves have been found, check collision of axis-aligned bounding boxes
      //        aabb_type bb1, bb2;
      //
      //        for (model_leaf_iterator it = leaf_data_.begin(); it != leaf_data_.end(); ++it) {
      //
      //          leaves_data_iterator lit = it->second.find(it1);
      //          if (lit != it->second.end())
      //            for (int i=0; i<lit->second.size(); ++i)
      //              bb1 += point_type(lit->second[i]);
      //
      //          lit = it->second.find(it2);
      //          if (lit != it->second.end())
      //            for (int i=0; i<lit->second.size(); ++i)
      //              bb2 += point_type(lit->second[i]);
      //        }
      //
      //        // assign velocities
      //        bb1.velocity_ = s1.velocity_;
      //        bb2.velocity_ = s2.velocity_;
      //
      //        // get time of intersection between aabbs
      //        Real ts = check_collision(bb1,bb2);
    }
    
    // traverse down the left hierarchy if its volume is bigger or if the
    // right volume is a leaf
    if ((!it1.is_leaf() && it1->measure() >= it2->measure()) || it2.is_leaf()) {
      
#ifdef DEBUG_MANAGER
      cout<<"    s1 >= s2 OR s2 is leaf"<<endl;
#endif
      iterator lit = it1.left();
      iterator rit = it1.right();
      assert(lit != null_);
      assert(rit != null_);
      check_collision(t, Dt, lit, it2, Bool2Type<Accelerate>());
      check_collision(t, Dt, rit, it2, Bool2Type<Accelerate>());
    }
    // traverse down the right hierarchy if its volume is bigger or if the
    // right volume is a leaf
    else if ((!it2.is_leaf() && it2->measure() > it1->measure()) || it1.is_leaf()) {
      
#ifdef DEBUG_MANAGER
      cout<<"    s1 < s2 OR s1 is leaf"<<endl;
#endif
      
      iterator lit = it2.left();
      iterator rit = it2.right();
      assert(lit != null_);
      assert(rit != null_);
      check_collision(t, Dt, it1, lit, Bool2Type<Accelerate>());
      check_collision(t, Dt, it1, rit, Bool2Type<Accelerate>());
    }
    else {
      // should never call this branch
      assert (false);
    }
  }
  
  
  
  // check collision considering acceleration
  // the equation to solve is
  //
  //   alpha t^4/4 + beta t^3 + gamma t^2 + delta t + epsilon = 0
  //
  // where alpha = (a.a)/4, beta = (a.v), gamma = a.s + v.v, delta = 2(s.v),
  // epsilon = s.s - r^2, a = a1-a2, v = v1-v2, s = c1-c2, r = r1-r2
  //
  template <class iterator>
  void check_collision(time_type t, time_type& Dt, iterator it1, iterator it2, Bool2Type<true>) {
    
    
    // check for collision
    if (*it1 & *it2)
      handle_collision(t, Dt, it1, it2);

    
    const point_type& v1 = it1->velocity_;
    const point_type& v2 = it2->velocity_;
    
    point_type a1 = it1->acceleration_;
    point_type a2 = it2->acceleration_;
    
    // TODO REMOVE HARDCODED ACCELERATION
    a1[0] = a1[1] = a2[0] = a2[1] = 0;
    a1[2] = a1[2] = -9.81;
    
    Real r = it1->radius() + it2->radius();
    point_type s = it2->center() - it1->center();
    point_type v = v2 - v1;
    point_type a = a2 - a1;
    
    // get coefficients
    Real alpha = (a*a)/4;
    Real beta = a*v;
    Real gamma = a*s + v*v;
    Real delta = 2*(s*v);
    Real epsilon = s*s - r*r;
    
    // obtain roots from quartic equation by calling the function such that
    // the coefficient for the quartic term is equal to one
    std::vector<Real> x(4, -1);
    
    cout<<"calling solve quartic with coefficients "<<endl;
    cout<<(beta/alpha)<<", "<<(gamma/alpha)<<", "<<(delta/alpha)<<", "<<(epsilon/alpha)<<endl;
    uint32_t roots = solve_quartic(beta/alpha, gamma/alpha, delta/alpha, epsilon/alpha,
                                   &x[0], &x[1], &x[2], &x[3]);
    
    // if there are roots, take the first one as an indication of the collision
    if (roots > 0) {
      
      Real tmin = std::numeric_limits<Real>::infinity();
      for (size_t i=0; i<x.size(); ++i)
        if (x[i] > 0)
          tmin = std::min(tmin, x[i]);
      if (tmin+t != std::numeric_limits<Real>::quiet_NaN()) {
        timer_.push(std::make_tuple(tmin+t, it1, it2));
        // change time step
        cout<<"APPROXIMATE COLLISION TIME -> t + tmin = "<<t<<" + "<<tmin<<" = "<<(tmin+t)<<endl;
        if (tmin == 0) {
          assert(false);
        }
        
        //        if (tmin < Dt) {
        //          cout<<"CHANGING TIME STEP TO "<<tmin<<endl;
        //          Dt = tmin;
        //        }
      }
    }
    
    cout<<" roots -> "<<roots<<endl;
    cout<<" values -> ";
    for (size_t i=0; i<x.size(); ++i)
      cout<<x[i]<<", ";
    cout<<endl;
  }
  
  // check collision for constant velocity
  template <class iterator>
  void check_collision(time_type t, time_type& Dt, iterator it1, iterator it2, Bool2Type<false> flag) {
    
    const volume_type& s1 = *it1;
    const volume_type& s2 = *it2;
    
    const point_type& v1 = s1.velocity_;
    const point_type& v2 = s2.velocity_;
    
    // vector between spheress
    point_type s = s2.center() - s1.center();
    // relative motion of s1 with respect to stationary s0
    point_type v = v2 - v1;
    // sum of radii
    value_type r = s1.radius() + s2.radius();
    value_type c = s*s - r*r;
    
#ifdef DEBUG_MANAGER
    cout<<"Checking collisiong between:"<<endl;
    cout<<"  "<<s1<<", velocity: "<<v1<<endl;
    cout<<"  "<<s2<<", velocity: "<<v2<<endl;
    cout<<"  Relative velocity: "<<v<<endl;
#endif
    
    // already intersecting, call recursively on branches
    if (c < 0.) {
      
#ifdef DEBUG_MANAGER
      cout<<"  Intersecting bounding volumes"<<endl;
#endif
      handle_collision(t, Dt, it1, it2);
      
    }
    
    value_type epsilon = 1e-8;
    value_type a = v*v;
    // if spheres not moving relative to each other
    if (a < epsilon) {
#ifdef DEBUG_MANAGER
      cout<<"    Objects not moving relative to each other"<<endl;
      cout<<"    a: "<<a<<endl;
#endif
      return;
    }
    
    value_type b = v*s;
    // if spheres not moving towards each other
    if (b >= 0.) {
#ifdef DEBUG_MANAGER
      cout<<"    Objects not moving towards each other"<<endl;
      cout<<"    b: "<<b<<endl;
#endif
      return;
    }
    
    value_type d = b*b - a*c;
    // if no real-valued root (d < 0), spheres do not intersect
    // otherwise add time to timer
    if (d >= 0.) {
      time_type ts = (-b - sqrt(d))/a;
      if (ts > 0.) {
#ifdef DEBUG_MANAGER
        cout<<"    Objects intersects at time "<<(t+ts)<<endl;
#endif
        timer_.push(std::make_tuple(t + ts, it1, it2));
      }
#ifdef DEBUG_MANAGER
      else {
        cout<<"    ts negative: "<<ts<<endl;
      }
#endif
      
    }
#ifdef DEBUG_MANAGER
    else {
      cout<<"    Objects do not intersect"<<endl;
      cout<<"    d: "<<d<<endl;
    }
#endif
  }
  
  bool check_penetration(leaves_data_iterator it1, leaves_data_iterator it2) {
    
    const coordinate_type& c1 = it1->second;
    const coordinate_type& c2 = it2->second;
    
    for (size_t i=0; i<c1.size(); ++i) {
      
      point_type r(c1[i]);
      
      // swith on the number of coordinates of the second container
      switch (c2.size()) {
          
          // triangle case
        case 3:
        {
          // form plane from three points
          
          // get element extreme points
          point_type o(c2[0]);
          point_type p(c2[1]);
          point_type q(c2[2]);
          
          Plane pi(o,p,q);
          
          // check for actual penetration using relative velocity
          if (check_penetration(r, it1->first->velocity_ - it2->first->velocity_, pi))
            return true;
        }
          break;
          
        default:
          cout<<"*** ERROR *** Function check_penetration in file "<<__FILE__<<", line "<<__LINE__;
          cout<<", is not implemented for node container of size "<<c2.size()<<endl;
          cout<<"*** ABORTING *** "<<endl;
          exit(1);
          break;
      }
    } // loop over coordinates of first container
    
    return false;
  }
  
  // intersects aabbs a and b moving with constant velocities.
  // On intersection, return time of first contact
  time_type check_collision(const aabb_type& a, const aabb_type& b) {
    
    // overlapping aabbs
    if (a & b) {
      cout<<"*** ERROR *** Colliding bounding boxes"<<endl;
      cout<<"              Aborting..."<<endl;
      exit(1);
    }
    
    // compute relative velocity
    const point_type v = b.velocity_ - a.velocity_;
    
    // initialize times of first and last contact
    time_type tf = 0., tl = 0.;
    
    // for each axis, determine the time of first and last contact
    for (int i=0; i<point_type::dim(); ++i) {
      
      if (v[i] < 0.) {
        if (b.max(i) < a.min(i))
          return inf; // non-intersecting
        if (a.max(i) < b.min(i))
          tf = std::max((a.max(i) - b.min(i)) / v[i], tf);
        if (b.max(i) > a.min(i))
          tl = std::min((a.min(i) - b.max(i)) / v[i], tl);
      }
      if (v[i] > 0.) {
        if (b.min(i) > a.max(i))
          return inf; // non-intersecting
        if (b.max(i) < a.min(i))
          tf = std::max((a.min(i) - b.max(i)) / v[i], tf);
        if (a.max(i) > b.min(i))
          tl = std::min((a.max(i) - b.min(i)) / v[i], tl);
      }
    }
    return tf;
  }
  
  /*! \brief Solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d = 0.
   *
   *  Solves the quartic equation. Returns the number of roots
   *  found. Roots are filled in ascending order.
   *  Code taken from the Ion Beam Simulator, which is distrubuted under
   *  the terms of the GNU General Public License as published by the Free
   *  Software Foundation; either version 2 of the License, or (at your
   *  option) any later version.
   */
  uint32_t solve_quartic(Real a, Real b, Real c, Real d,
                         Real *x0, Real *x1, Real *x2, Real *x3 );
  
  
  friend std::ostream& operator<<(std::ostream& os, const Contact_model_manager2& mm) {
    
    os<<"Contact model manager info:"<<endl;
    os<<"  models: "<<mm.models_.size()<<endl;
    size_t i=0;
    const_forest_iterator tit = mm.forest_.begin();
    for (const_model_iterator it = mm.models_.begin(); it != mm.models_.end(); ++it) {
      os<<"\tmodel "<<++i<<" memory address: "<<*it<<endl;
      os<<"\tmodel: "<<**it<<endl;
      os<<"\ttree: ";
      print_mathematica(**tit++);
    }
    return os;
  }
};

__END_AKANTU__

#endif /* __AKANTU_MODEL_MANAGER_HH__ */
