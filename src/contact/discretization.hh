/**
 * @file   discretization.hh
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @date   Mon Jan 23 09:00:00 2012
 *
 * @brief  contact discretization classes
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

#ifndef __AKANTU_DISCRETIZATION_HH__
#define __AKANTU_DISCRETIZATION_HH__

#include <cassert>
#include <stdio.h>
#include <stdlib.h>

#ifdef AKANTU_USE_IOHELPER
#include <io_helper.hh>
#endif

#include "aka_common.hh"
#include "aka_types.hh"
#include "contact_common.hh"
#include "solid_mechanics_model.hh"
#include "friction.hh"

using std::endl;
using std::cout;

__BEGIN_AKANTU__



class Discretization_base {
  //! Support the visitor design pattern
  virtual void accept(Contact_discretization_visitor&) = 0;
};


template <int>
class Contact_discretization;

template <>
class Contact_discretization<N2N_t> : public Discretization_base {
  
  template <int> friend class Contact_scheme;
  
private:
  
  // pair structure
  struct N2N_pair {
    
    typedef UInt node_id;
    
    node_id master_, slave_;
    bool contact_, stick_;
    vector_type n_;
    
    N2N_pair(node_id m, node_id s)
    : master_(m), slave_(s), contact_(true), stick_(true), n_() {}
    
    N2N_pair(node_id m, node_id s, bool c, bool st) 
    : master_(m), slave_(s), contact_(c), stick_(st), n_() {}
    
    N2N_pair(std::istream& is) {
      is>>master_;
      is>>slave_;
      is>>contact_;
      is>>stick_;
      // NOTE: normal is not read from a file
    }
  };
  
public:
  
  /* iterators and pairing list */
  typedef std::list<N2N_pair> pairing_list;
  typedef pairing_list::iterator pairing_iterator;
  typedef pairing_list::const_iterator pairing_const_iterator;
  
//private:
  
  pairing_list node_pairs_;
  
public:
  
  //! Default constructor
  Contact_discretization() : node_pairs_() {}
  
  //! Node list parameter constructor
  template <class node_list>
  Contact_discretization(const node_list& masters, 
                         const node_list& slaves) : node_pairs_()
  { add_pairs(masters, slaves); }
  
  //! Construct the object from a file
  Contact_discretization(const std::string& str) {
    std::ifstream ifs(str.c_str());
    if (!ifs) {
      std::string err("Could not create node-to-node discretization from file " + str);
      AKANTU_EXCEPTION(err);
    }
    ifs>>*this;
  }
  
  void set_normals(const vector_type& n) {
    for (pairing_iterator it = node_pairs_.begin();
         it != node_pairs_.end(); ++it)
      it->n_ = n;
  }
  
  std::pair<UInt,UInt> get_contact_info() const {
    
    UInt contact = 0, stick = 0;
    for (pairing_const_iterator it = node_pairs_.begin(); it != node_pairs_.end(); ++it)
      if (it->contact_) {
        ++contact;
        if (it->stick_)
          ++stick;
      }
    return std::pair<UInt,UInt>(contact, stick);
  }
  
  //! Member function to support the visitor design pattern.
  void accept(Contact_discretization_visitor& guest)
  { guest.GenericVisit(*this); }
  
  
  //! Add node-to-node pairs
  template <class node_list>
  void add_pairs(const node_list& masters, const node_list& slaves) {
    
    if (masters.size() != slaves.size())
      AKANTU_EXCEPTION("In Discretization_N2N constructor, "
                       << "master and slave containers of different size");
    
    typedef typename node_list::const_iterator node_iterator;
    for (node_iterator it1 = masters.begin(), it2 = slaves.begin();
         it1 != masters.end(); ++it1, ++it2)
      node_pairs_.push_back(N2N_pair(*it1, *it2));
  }
  
  friend std::ostream& operator<<(std::ostream& os, const Contact_discretization& d) {
    
    os<<"# Node-to-node discretization pair information:"<<endl;
    for (pairing_const_iterator it = d.node_pairs_.begin(); it != d.node_pairs_.end(); ++it)
      os<<it->master_<<" "<<it->slave_<<" "<<it->contact_<<" "<<it->stick_<<endl;
    return os;
  }
  
  friend std::istream& operator>>(std::istream& is, Contact_discretization& d) {
    
    while(is.good()) {
      
      std::string line;
      getline(is, line);
      
      // look for header
      if (line.find("Node-to-node") != std::string::npos)
        // read pairs
        while (is.good())
          d.node_pairs_.push_back(N2N_pair(is));
    }
    return is;
  }
  
};






template <Friction_type friction>
class Contact_MPC : public Contact_friction<friction> {
  
  typedef Contact_friction<friction> friction_type;
  
  /* Node to Node structure */
  struct N2N {
    
    UInt master, slave;
    bool contact, stick;
    
    N2N(UInt m, UInt s)
    : master(m), slave(s), contact(true), stick(true) {}
    
    N2N(UInt m, UInt s, bool c, bool st) 
    : master(m), slave(s), contact(c), stick(st) {}
    
  };
  
  /* iterators and pairing list */
  typedef std::list<N2N> pairing_list;
  typedef typename pairing_list::iterator pairing_iterator;
  typedef typename pairing_list::const_iterator pairing_const_iterator;
  
  /* members */
  const SolidMechanicsModel & model;
  pairing_list node_pairs;
  
  /* functions */
public:
  
  template <class node_set>
  Contact_MPC(const SolidMechanicsModel & model, 
              const node_set& masters, 
              const node_set& slaves);
  
  Contact_MPC(const SolidMechanicsModel & model);
  
  UInt getNbContactNodes() const;
  UInt getNbStickNodes() const;
  
  void frictionPred();
  void syncAcceleration();
  
  void dumpRestart(const std::string&);
  void readRestart(const std::string&);
  
  friend std::ostream& operator<<(std::ostream& os, const Contact_MPC& c) {
    os<<"Contact multiple point constraint node info:"<<endl;
    for (typename Contact_MPC::pairing_const_iterator it = c.node_pairs.begin(); it != c.node_pairs.end(); ++it) {
      os<<"  ("<<it->master<<"-"<<it->slave<<"), contact: "<<it->contact<<", stick: "<<it->stick<<endl;
    }
    return os;
  }
};


template <Friction_type friction>
template <class node_set>
Contact_MPC<friction>::Contact_MPC(const SolidMechanicsModel & model, 
                                   const node_set& masters, 
                                   const node_set& slaves) : model(model), node_pairs() {
  
  assert(masters.size() == slaves.size());
  
  typedef typename node_set::const_iterator node_iterator;
  for (node_iterator it1 = masters.begin(), it2 = slaves.begin();
       it1 != masters.end(); ++it1, ++it2) {
    
    node_pairs.push_back(N2N(*it1, *it2));
    
  }
}


template <Friction_type friction>
Contact_MPC<friction>::Contact_MPC(const SolidMechanicsModel & model) : model(model), node_pairs() {}

template <Friction_type friction>
UInt Contact_MPC<friction>::getNbContactNodes() const {
  UInt nb_contact = 0;
  
  for (typename Contact_MPC::pairing_const_iterator it = node_pairs.begin(); 
       it != node_pairs.end(); ++it) {
    if (it->contact) nb_contact++;
  }
  return nb_contact;
}

template <Friction_type friction>
UInt Contact_MPC<friction>::getNbStickNodes() const {
  UInt nb_sticks = 0;
  
  for (typename Contact_MPC::pairing_const_iterator it = node_pairs.begin(); 
       it != node_pairs.end(); ++it) {
    if (it->contact && it->stick) ++nb_sticks;
  }
  return nb_sticks;
}

template <Friction_type friction>
void Contact_MPC<friction>::frictionPred() {
  const Vector<Real> & residual = model.getResidual();
  
  for (typename Contact_MPC::pairing_iterator it = node_pairs.begin(); 
       it != node_pairs.end(); ++it) {
    if (it->contact && it->stick) {
      
      Real normal_pressure = 0.5 * (fabs(residual(it->slave,1)) + fabs(residual(it->master,1)));
      Real rel_tang_traction = fabs(residual(it->slave,0) - residual(it->master,0));
      
      if (rel_tang_traction > this->computeForceMagnitude(normal_pressure))
        it->stick = false;
      else
        it->stick = true;
    }
    // else if slip
    else if (it->contact && !it->stick) {
      
      
      
    }
  }
}

template <Friction_type friction>
void Contact_MPC<friction>::syncAcceleration() {
  Vector<Real> & acceleration = model.getIncrementAcceleration();
  const Vector<Real> & mass = model.getMass();
  const Vector<Real> & residual = model.getResidual();
  
  for (typename Contact_MPC::pairing_const_iterator it = node_pairs.begin(); 
       it != node_pairs.end(); ++it) {
    if (it->contact) {
      Real common_n_acc = (residual(it->slave,1) + residual(it->master,1))/(mass(it->slave,1) + mass(it->master,1));
      //	cout << common_n_acc << " ";
      acceleration(it->slave,1) = common_n_acc;
      acceleration(it->master,1) = common_n_acc;
      
      if (it->stick) {
        Real common_t_acc = (residual(it->slave,0) + residual(it->master,0))/(mass(it->slave,0) + mass(it->master,0));
        acceleration(it->slave,0) = common_t_acc;
        acceleration(it->master,0) = common_t_acc;
      }
    }
  }
}

template <Friction_type friction>   
void Contact_MPC<friction>::dumpRestart(const std::string& filename) {
  
  std::ofstream out_restart;
  out_restart.open(filename.c_str());
  out_restart << "%master slave contact stick" << std::endl;
  
  for (typename Contact_MPC::pairing_const_iterator it = node_pairs.begin(); 
       it != node_pairs.end(); ++it) {
    out_restart << it->master << " " << it->slave << " "<< it->contact << " " << it->stick << endl;
  }
  
  out_restart.close();
}

template <Friction_type friction>
void Contact_MPC<friction>::readRestart(const std::string& filename) {
  
  std::ifstream infile(filename.c_str());
  
  if (!infile)
    cout<<"*** ERROR *** File "<<filename<<" could not be opened. In file "<<__FILE__<<", line "<<__LINE__<<endl;
  
  std::string line;
  
  while(infile.good()) {
    char c = infile.peek();
    if (c == '#')
      continue;
    
    getline(infile, line);
    UInt master, slave;
    bool contact, stick;
    std::stringstream sstr(line);
    sstr>>master; sstr>>slave; sstr>>contact; sstr>>stick;
    
    node_pairs.push_back(N2N(master, slave, contact, stick));
  }
}

__END_AKANTU__

#endif /* __AKANTU_DISCRETIZATION_HH__ */
