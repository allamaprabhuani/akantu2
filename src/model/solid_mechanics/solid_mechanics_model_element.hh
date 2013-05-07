/**
 * @file   solid_mechanics_model_element.hh
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Aug 22 12:14:00 2012
 *
 * @brief  elements for solid mechanics models
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


#ifndef AKANTU_SOLID_MECHANICS_MODEL_ELEMENT_HH
#define AKANTU_SOLID_MECHANICS_MODEL_ELEMENT_HH

#include <array/expr.hpp>

#include "mesh.hh"
#include "aka_error.hh"
#include "solid_mechanics_model.hh"


__BEGIN_AKANTU__


typedef array::vector_type<Real> vector_type;
typedef array::matrix_type<Real> matrix_type;




template <>
class ModelElement<SolidMechanicsModel> : public Element {
  
  typedef Element element_base;

  UInt *connectivity_;               //!< Ponter to connectivity array
  SolidMechanicsModel &model_;       //!< Reference to model
  ElementType type_;                 //!< Element type
  
public:
  
  typedef SolidMechanicsModel model_type;
  
  ModelElement(SolidMechanicsModel& m, ElementType type, UInt id, GhostType gt = _not_ghost)
  : element_base(type, id, gt), model_(m), type_(type) {
    connectivity_ = &m.getMesh().getConnectivity(type, gt)(id);
  }
  
  model_type& model() { return model_; }
  
  UInt& node(UInt n) {
    AKANTU_DEBUG_ASSERT(n < Mesh::getNbNodesPerElement(type_),
                        "Node "<<n<<" is larger than element number of nodes: "<<Mesh::getNbNodesPerElement(type_));
    return connectivity_[n];
  }
  
  UInt node(UInt n) const {
    AKANTU_DEBUG_ASSERT(n < Mesh::getNbNodesPerElement(type_),
                        "Node "<<n<<" is larger than element number of nodes: "<<Mesh::getNbNodesPerElement(type_));
    return connectivity_[n];
  }
  
  
  // vector of pointers to nodes' first coordinates
  std::vector<const Real*> coordinates() {
    UInt nb_nodes = Mesh::getNbNodesPerElement(this->type);
    const Array<Real> &position = model_.getCurrentPosition();
    std::vector<const Real*> coord(nb_nodes);
    for (size_t i=0; i<nb_nodes; ++i)
      coord[i] = &position(connectivity_[i]);
    return coord;
  }
  
  // barycenter
  vector_type barycenter() const {
    
    typedef typename vector_type::value_type value_type;
    
    UInt nb_nodes = Mesh::getNbNodesPerElement(this->type);
    const Array<Real> &position = model_.getCurrentPosition();
    vector_type sum(model_.getSpatialDimension());
    for (size_t i=0; i<nb_nodes; ++i) {
      Real * p = const_cast<Real*>(&position(connectivity_[i]));
      sum += vector_type(model_.getSpatialDimension(), p);
    }

    return (1./static_cast<value_type>(nb_nodes)) * sum;
  }
  
  
  // mass
  vector_type getMass(UInt nid) {
    UInt n = node(nid);
    Array<Real> &mass = model_.getMass();
    return vector_type(&mass(n), model_.getSpatialDimension());
  }
  
  // mass for const objects
  const vector_type getMass(UInt nid) const {
    UInt n = node(nid);
    Array<Real> &mass = model_.getMass();
    return vector_type(&mass(n), model_.getSpatialDimension());
  }
  
  // displacement
  vector_type getDisplacement(UInt nid) {
    UInt n = node(nid);
    Array<Real> &displacement = model_.getDisplacement();
    return vector_type(&displacement(n), model_.getSpatialDimension());
  }
  
  // displacement for const objects
  const vector_type getDisplacement(UInt nid) const {
    UInt n = node(nid);
    Array<Real> &displacement = model_.getDisplacement();
    return vector_type(&displacement(n), model_.getSpatialDimension());
  }
  
  // velocity
  vector_type getVelocity(UInt nid) {
    UInt n = node(nid);
    Array<Real> &velocity = model_.getVelocity();
    return vector_type(&velocity(n), model_.getSpatialDimension());
  }
  
  // velocity for const objects
  const vector_type getVelocity(UInt nid) const {
    UInt n = node(nid);
    Array<Real> &velocity = model_.getVelocity();
    return vector_type(&velocity(n), model_.getSpatialDimension());
  }
  
  // acceleration
  vector_type getAcceleration(UInt nid) {
    UInt n = node(nid);
    Array<Real> &acceleration = model_.getAcceleration();
    return vector_type(&acceleration(n), model_.getSpatialDimension());
  }
  
  // acceleration for const objects
  const vector_type getAcceleration(UInt nid) const {
    UInt n = node(nid);
    Array<Real> &acceleration = model_.getAcceleration();
    return vector_type(&acceleration(n), model_.getSpatialDimension());
  }
  
  // position (location + displacement)
  vector_type getCurrentPosition(UInt nid) {
    UInt n = node(nid);
    const Array<Real> &position = model_.getCurrentPosition();
    Real * p = const_cast<Real*>(&position(n));
    return vector_type(model_.getSpatialDimension(), p);
  }
  
  // position (location + displacement) for const objects
  const vector_type getCurrentPosition(UInt nid) const {
    UInt n = node(nid);
    const Array<Real> &position = model_.getCurrentPosition();
    Real * p = const_cast<Real*>(&position(n));
    return vector_type(model_.getSpatialDimension(), p);
  }
  
  // residual
  vector_type getResidual(UInt nid) {
    UInt n = node(nid);
    const Array<Real> & residual = model_.getResidual();
    Real * p = const_cast<Real*>(&residual(n));
    return vector_type(model_.getSpatialDimension(), p);
  }

  // residual for const objects
  const vector_type getResidual(UInt nid) const {
    UInt n = node(nid);
    const Array<Real> & residual = model_.getResidual();
    Real * p = const_cast<Real*>(&residual(n));
    return vector_type(model_.getSpatialDimension(), p);
  }
  
  // momentum
  vector_type getMomentum(UInt nid) {
    UInt n = node(nid);
    Array<Real> &velocity = model_.getVelocity();
    Array<Real> &mass = model_.getMass();
    vector_type p(model_.getSpatialDimension());
    for (size_t i=0; i<p.size(); ++i)
      p[i] = mass(n,i)*velocity(n,i);
    return p;
  }
  
  // momentum for const objects
  const vector_type getMomentum(UInt nid) const {
    return const_cast<ModelElement&>(*this).getMomentum(nid);
  }

};


__END_AKANTU__


#endif /* AKANTU_SOLID_MECHANICS_MODEL_ELEMENT_HH */
