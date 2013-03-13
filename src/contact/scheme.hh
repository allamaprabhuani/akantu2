/**
 * @file   scheme.hh
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @date   Mon Jan 23 09:00:00 2012
 *
 * @brief  contact scheme classes
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

#ifndef __AKANTU_SCHEME_HH__
#define __AKANTU_SCHEME_HH__

#include "aka_common.hh"
#include "contact_common.hh"
#include "discretization.hh"
#include "friction.hh"

__BEGIN_AKANTU__



template <int>
class Contact_scheme;

//! Explicit scheme
template<>
class Contact_scheme<Explicit_t> : public Contact_discretization_visitor {
  
  
  //    // Visit discretizations
  //    void Visit(Contact_discretization<N2N_t>& cd) {
  //        
  //        cout<<"VISITING NODE TO NODE DISCRETIZATION1!!!!"<<endl;
  //    }
  
public:

  // WORKING WITHOUT FRICTION
//  template <class model_type>
//  void frictionPred(const model_type& model,
//                    Contact_discretization<N2N_t>& disc,
//                    Contact_friction<Coulomb_t>& friction) {
//    
//    typedef Contact_discretization<N2N_t> discretization_type;
//    
//    const Array<Real> &R = model.getResidual(); // residual
//    
//    UInt d = model.getSpatialDimension();
//    
//    matrix_type I = eye(d);
//    
//    for (typename discretization_type::pairing_iterator it = disc.node_pairs_.begin(); 
//         it != disc.node_pairs_.end(); ++it) {
//      
//      if (it->contact_ && it->stick_) {
//        
//        const vector_type& nm = it->n_;
//        vector_type ns = -1.*nm;
//        
//        vector_type fm(d), fs(d);
//        for (UInt i=0; i<d; ++i) {
//          fm[i] = R(it->master_,i);
//          fs[i] = R(it->slave_,i);
//        }
//
//        // get normal pressure 
//        Real np1 = transpose(fm)*nm;
//        Real np2 = transpose(fs)*ns;
//        Real normal_pressure = 0.5 * std::abs(np1) + std::abs(np2); // n_m = -n_s
//        
//        // TODO remove in contact if negative
//        
//        // obtain projection operators
//        matrix_type Pm = I - nm * transpose(nm);
//        matrix_type Ps = I - ns * transpose(ns);
//        
//        // get tangential vector
//        vector_type tau = Pm*fm - Ps*fs;
//        Real rel_tang_traction = tau.norm();
//        
//        // set stick to false if tangential traction is greater than friction force
//        if (rel_tang_traction > friction.computeForceMagnitude(normal_pressure))
//          it->stick_ = false;
//
//      }
//      // else if slip
//      else if (it->contact_ && !it->stick_) {
//        
//        // TODO write slip case
//        
//      }
//    }
//  }
  
  
  template <class iterator>
  void friction_against_residual(iterator it,
                                 const vector_type& rm,
                                 const vector_type& rs,
                                 Real fric_strength,
                                 Real* residual) {
    UInt d = rm.size();
    matrix_type I = eye(d);

    const vector_type& nm = it->n_;
    vector_type ns = -1.*nm;

    // obtain projection operators
    matrix_type Pm = I - nm * transpose(nm);
    matrix_type Ps = I - ns * transpose(ns);
    
    // get tangential vectors
    vector_type rm_tau = Pm*rm;
    vector_type rs_tau = Ps*rs;
    Real rm_tau_norm = rm_tau.norm();
    Real rs_tau_norm = rs_tau.norm();
    
    // set stick to false if tangential traction is greater than friction force
    if (rm_tau_norm > fric_strength || rs_tau_norm > fric_strength) {
      // set slip condition
      it->stick_ = false;
      // normalize tangential vectors
      rm_tau /= rm_tau_norm;
      rs_tau /= rs_tau_norm;
      // add frictional forces in the opposite direction
      
      // add frictional force to residual due to sliding
      for (UInt i=0; i<d; ++i) {
        residual[it->master_*d + i] -= fric_strength * rm_tau(i);
        residual[it->slave_*d + i] -= fric_strength * rs_tau(i);
      }
    } else {
      it->stick_ = true;
    // TODO check accelerations to make sure they are accurate after they stick
    }
    
  }
  
  
  template <class model_type>
  void frictionPred(const model_type& model,
                    Contact_discretization<N2N_t>& disc,
                    Contact_friction<Coulomb_t>& friction) {
    
    typedef Contact_discretization<N2N_t> discretization_type;
    
    const Array<Real> &R = model.getResidual(); // residual
    Real *residual = R.storage();
    
    UInt d = model.getSpatialDimension();
    Real tol = Math::getTolerance();

    matrix_type I = eye(d);
    
    for (typename discretization_type::pairing_iterator it = disc.node_pairs_.begin(); 
         it != disc.node_pairs_.end(); ++it) {
      
      if (!it->contact_) {
        // set to default stick state in case nodes enter in contact again
        it->stick_ = true;
        continue;
      }
      
      const vector_type& nm = it->n_;
      vector_type ns = -1.*nm;
      
      vector_type rm(d), rs(d);
      for (UInt i=0; i<d; ++i) {
        rm[i] = R(it->master_,i);
        rs[i] = R(it->slave_,i);
      }
      
      // get normal pressure 
      Real np1 = transpose(rm)*nm;
      Real np2 = transpose(rs)*ns;

      // obtain normal pressure
      Real normal_pressure = 0.5 * (np1 + np2);
      
      // compute frictional strength for the level of normal pressure
      Real fric_strength = friction.computeForceMagnitude(normal_pressure);
      
      // sync acceleration if in contact and stick
      if (it->stick_) {
        
        // remove contact if negative contact pressure
        if (np1 + np2 < -tol) {
          it->contact_ = false;
          cout<<"*** INFO *** removing contact"<<endl;
          // do not do anything to contact forces, they are zero
          continue;
        } else {
          
          // add friction against residual
          friction_against_residual(it, rm, rs, fric_strength, residual); 
        }
      }
      // else if slip
      else {
        
        // compute velocities
        const Array<Real> &V = model.getVelocity();
        
        vector_type vm(d), vs(d);
        for (UInt i=0; i<d; ++i) {
          vm[i] = V(it->master_,i);
          vs[i] = V(it->slave_,i);
        }
        
        // obtain relative velocity vector
        vector_type v_rel = vm-vs;
        
        // if there is a non-zero relative velocity
        if (v_rel.norm() > tol) {
          
          // normalize velocities to obtain unit vector
          vm /= vm.norm();
          vs /= vs.norm();
          
          // add frictional force agains the direction of the velocity vector
          for (UInt i=0; i<d; ++i) {
            residual[it->master_*d + i] -= fric_strength * vm(i);
            residual[it->slave_*d + i] -= fric_strength * vm(i);
          }
        } else {
          
          // add friction agains residual
          friction_against_residual(it, rm, rs, fric_strength, residual);
          
          cout<<"*** INFO *** No velocity after sliding"<<endl;


        }
      }
    }
  }
  
  
  // WORKING WITHOUT FRICTION
//  template <class model_type>
//  void syncAcceleration(const model_type& model,
//                        Contact_discretization<N2N_t>& disc,
//                        Contact_friction<Coulomb_t>& friction) {
//    
//    typedef Contact_discretization<N2N_t> discretization_type;
//    
//    Array<Real> &A = model.getIncrementAcceleration();
//    const Array<Real> &M = model.getMass();
//    const Array<Real> &R = model.getResidual();
//    
//    // problem dimension
//    UInt d = model.getSpatialDimension();
//    
//    // identity matrix
//    matrix_type I = eye(d);
//    
//    // loop over interface pairs
//    for (typename discretization_type::pairing_iterator it = disc.node_pairs_.begin(); 
//         it != disc.node_pairs_.end(); ++it) {
//      
//      // synchronize only if in contact
//      if (it->contact_) {
//        
//        // if in contact, a portion of the acceleration will be synchronized
//        // between the nodes (only the normal component if no-stick condition)
//        // NOTE: it is assumed that the mass matrix contains already the mass
//        // in the global coordinates, otherwise a rotation matrix considering
//        // the normal has to be used
//        
//        // first consider both nodes as a single node
//        vector_type r(d);
//        matrix_type M_inv(d,d);
//        for (UInt i=0; i<d; ++i) {
//          r[i] = R(it->master_,i) + R(it->slave_,i);
//          M_inv[i][i] = 1./(M(it->master_,i) + M(it->slave_,i));
//        }
//        
//        // obtain global acceleration vector
//        vector_type a = M_inv*r;
//        
//        // if stick, synchorinize all acceleration components
//        if (it->stick_) {
//          
//          // synchronize acceleration
//          for (UInt i=0; i<d; ++i)
//            A(it->slave_,i) = A(it->master_,i) = a(i);
//          
//          // synchronize only normal projection
//        } else {
//          
//          const vector_type& nm = it->n_;
//          vector_type ns = -1.*nm;
//          
//          // obtain normal component
//          Real s = transpose(a)*nm;
//          vector_type a_n = s*nm;
//          
//          // consider nodes individually to obtain tangential components
//          vector_type rm(d), rs(d);
//          matrix_type Mm_inv(d,d), Ms_inv(d,d);
//          for (UInt i=0; i<d; ++i) {
//            rm[i] = R(it->master_,i);
//            rs[i] = R(it->slave_,i);
//            Mm_inv[i][i] = 1./M(it->master_,i);
//            Ms_inv[i][i] = 1./M(it->slave_,i);
//          }
//          
//          // obtain projection operators
//          matrix_type Pm = I - nm * transpose(nm);
//          matrix_type Ps = I - ns * transpose(ns);
//          
//          // get acceleration vectors
//          vector_type am = Mm_inv*rm;
//          vector_type as = Ms_inv*rs;
//          
//          // obtain tangential accelerations
//          vector_type am_tau = Pm*am;
//          vector_type as_tau = Ps*as;
//          
//          // synchronize acceleration
//          for (UInt i=0; i<d; ++i) {
//            A(it->master_,i) = am_tau(i) + a_n(i);
//            A(it->slave_,i) = as_tau(i) + a_n(i);
//          }
//          
//        } // no stick
//      } // if contact
//    } // loop over interface pairs
//  }
  
  
  template <class model_type>
  void syncAcceleration(const model_type& model,
                        Contact_discretization<N2N_t>& disc,
                        Contact_friction<Coulomb_t>& friction) {
    
    typedef Contact_discretization<N2N_t> discretization_type;
    
    Array<Real> &A = model.getIncrementAcceleration();
    const Array<Real> &M = model.getMass();
    const Array<Real> &R = model.getResidual();
    
    // problem dimension
    UInt d = model.getSpatialDimension();
    
    // identity matrix
    matrix_type I = eye(d);
    
    // loop over interface pairs
    for (typename discretization_type::pairing_iterator it = disc.node_pairs_.begin(); 
         it != disc.node_pairs_.end(); ++it) {
      
      // synchronize only if in contact
      if (it->contact_) {
        
        // if in contact, a portion of the acceleration will be synchronized
        // between the nodes (only the normal component if no-stick condition)
        // NOTE: it is assumed that the mass matrix contains already the mass
        // in the global coordinates, otherwise a rotation matrix considering
        // the normal has to be used
        
        // first consider both nodes as a single node
        vector_type r(d);
        matrix_type M_inv(d,d);
        for (UInt i=0; i<d; ++i) {
          r[i] = R(it->master_,i) + R(it->slave_,i);
          M_inv[i][i] = 1./(M(it->master_,i) + M(it->slave_,i));
        }
        
        // obtain global acceleration vector
        vector_type a = M_inv*r;
        
        // if stick, synchorinize all acceleration components
        if (it->stick_) {
          
          // synchronize acceleration
          for (UInt i=0; i<d; ++i)
            A(it->slave_,i) = A(it->master_,i) = a(i);
          
          // synchronize only normal projection
        } else {
          
          const vector_type& nm = it->n_;
          vector_type ns = -1.*nm;
          
          // obtain normal component
          Real s = transpose(a)*nm;
          vector_type a_n = s*nm;
          
          // consider nodes individually to obtain tangential components
          vector_type rm(d), rs(d);
          matrix_type Mm_inv(d,d), Ms_inv(d,d);
          for (UInt i=0; i<d; ++i) {
            rm[i] = R(it->master_,i);
            rs[i] = R(it->slave_,i);
            Mm_inv[i][i] = 1./M(it->master_,i);
            Ms_inv[i][i] = 1./M(it->slave_,i);
          }
          
          // obtain projection operators
          matrix_type Pm = I - nm * transpose(nm);
          matrix_type Ps = I - ns * transpose(ns);
          
          // get acceleration vectors
          vector_type am = Mm_inv*rm;
          vector_type as = Ms_inv*rs;
          
          // obtain tangential accelerations
          vector_type am_tau = Pm*am;
          vector_type as_tau = Ps*as;
          
          // synchronize acceleration
          for (UInt i=0; i<d; ++i) {
            A(it->master_,i) = am_tau(i) + a_n(i);
            A(it->slave_,i) = as_tau(i) + a_n(i);
          }
          
        } // no stick
      } // if contact
    } // loop over interface pairs
  }

  
};


//! Implicit scheme
template<>
class Contact_scheme<Implicit_t> {
  
  
  
};



__END_AKANTU__

#endif /* __AKANTU_SCHEME_HH__ */
