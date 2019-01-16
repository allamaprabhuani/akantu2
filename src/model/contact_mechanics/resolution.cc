/**
 * @file   resolution.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Jan 7 2019
 * @date last modification: Mon Jan 7 2019
 *
 * @brief  Implementation of common part of the contact resolution class
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "resolution.hh"
#include "contact_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
Resolution::Resolution(ContactMechanicsModel & model, const ID & id)
  : Memory(id, model.getMemoryID()), Parsable(ParserType::_contact_resolution, id),
    fem(model.getFEEngine()),
    name(""), model(model),
    spatial_dimension(model.getMesh().getSpatialDimension()){

  AKANTU_DEBUG_IN();
  
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Resolution::~Resolution() = default;

/* -------------------------------------------------------------------------- */
void Resolution::initialize() {
  registerParam("name", name, std::string(), _pat_parsable | _pat_readable); 
  registerParam("mu", mu, Real(0.), _pat_parsable | _pat_modifiable,
		"Friciton Coefficient");
}

/* -------------------------------------------------------------------------- */
void Resolution::printself(std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  std::string type = getID().substr(getID().find_last_of(':') + 1);

  stream << space << "Contact Resolution " << type << " [" << std::endl;
  Parsable::printself(stream, indent);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void Resolution::assembleInternalForces(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const Array<Int> & equation_array =
    model.getDOFManager().getEquationsNumbers();
 
  auto & internal_force =
    const_cast<Array<Real> &>(model.getInternalForce());

  auto & contact_map = model.getContactMap();
  
  const auto slave_nodes = model.getMesh().getNodeGroup(name);
  
  for (auto & slave: slave_nodes) {

    auto & master     = contact_map[slave].master;
    auto & gap        = contact_map[slave].gap;
    auto & projection = contact_map[slave].projection;

    UInt nb_nodes_master = Mesh::getNbNodesPerElement(master.type)
    
    Vector<Real> shapes(nb_nodes_master);
    fem.computeShapes(projection, master, master.type, shapes, ghost_type);
    
    Matrix<Real> shapes_derivatives(spatial_dimension - 1, nb_nodes_master);
    fem.computeShapeDerivatives(projection, master, master.type, shapes_derivatives, ghost_type);

    const auto & connectivity = contact_map[slave].connectivity;
    Vector<Real> elementary_force(connectivity.size() * spatial_dimension);

    Array<Real> * tangents =
      new Array<Real>(spatial_dimension - 1, spatial_dimension, "surface_tangents");

    Array<Real> * global_coords =
      new Array<Real>(nb_nodes_master, spatial_dimension);

    computeCoordinates(master.type, *global_coords);
    computeTangents(shapes_derivatives, *global_coords, *tangents);

    Matrix<Real> surface_matrix(spatial_dimension - 1, spatial_dimension - 1);
    computeSurfaceMatrix(*tangents, surface_matrix);

    computeN(*n, shapes, normal);
    computeNormalForce(elementary_force, *n, gap);

    computeTalpha(*t_alpha,  shapes,            *tangents);
    computeNalpha(*n_alpha, *shapes_derivatives, normal);
    computeDalpha(*d_alpha, *n_alpha, *t_alpha, surface_matrix);
    computeFrictionForce(elementary_force, *d_alpha, gap);
         
    for (UInt i = 0; i < connectivity.size(); ++i) {
      for (UInt j = 0; j < spatial_dimension; ++j) {	
	UInt offset_node = connectivity(i) * spatial_dimension + j;
	auto & equation_num = equation_array(offset_node);
	internal_force(equation_num) += elementary_force(i + j);
      }
    }    
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Resolution::assembleStiffnessMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const auto slave_nodes = model.getMesh().getNodeGroup(name);
  auto & contact_map = model.getContactMap();
  
  for (auto & slave: slave_nodes) {

    auto & master     = contact_map[slave].master;
    auto & gap        = contact_map[slave].gap;
    auto & projection = contact_map[slave].projection;
    
    Vector<Real> shapes(master.nb_nodes);
    fem.computeShapes(projection, master, master.type, shapes, ghost_type);

    Vector<Real> shapes_derivatives(master.nb_nodes * spatial_dimension);
    fem.computeShapeDerivatives(projection, master, master.type, shapes_derivatives, ghost_type);

    const auto & connectivity = contact_map[slave].connectivity;
    Matrix<Real> elementary_stiffness(connectivity.size() * spatial_dimension,
				      connectivity.size() * spatial_dimension);

    Array<Real> * tangents =
      new Array<Real>(spatial_dimension - 1, spatial_dimension, "surface_tangents");

    Array<Real> * global_coords =
      new Array<Real>(nb_nodes_master, spatial_dimension);

    computeCoordinates(master.type, *global_coords);
    computeTangents(shapes_derivatives, *global_coords, *tangents);

    Matrix<Real> surface_matrix(spatial_dimension - 1, spatial_dimension - 1);
    computeSurfaceMatrix(*tangents, surface_matrix);

    computeN(     *n,        shapes,             normal);
    computeTalpha(*t_alpha,  shapes,            *tangents);
    computeNalpha(*n_alpha, *shapes_derivatives, normal);
    computeDalpha(*d_alpha, *n_alpha, *t_alpha,  surface_matrix);

    computeTangentModuli(*n, *n_alpha, *t_alpha, *d_alpha, gap);
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Resolution::computeTangents(Matrix<Real> & shapes_derivatives, Array<Real> & global_coords,
				 Array<Real> & tangents) {

  UInt i = 0;
  for (auto && values : zip(make_view(tangents, spatial_dimension))) {
    auto & tangent = std::get<0>(values);
    for (UInt n : arange(global_coords.nb_components)) {
      tangent += shapes_derivaties(n, i) * global_coords(n);
    }
    ++i;
  }

}

/* -------------------------------------------------------------------------- */
void Resolution::computeSurfaceMatrix(Array<Real> & tangents, Matrix<Real> & surface_matrix) {

  for (UInt i : arange(spatial_dimension - 1)) {
    for (UInt j : arange(spatial_dimension -1 )) {
      surface_matrix(i, j) = tangents(i) * tangents(j);
    }
  }

  inverse(surface_matrix);
}

/* -------------------------------------------------------------------------- */
void Resolution::computeN(Array<Real> & n, Vector<Real> & shapes, Vector<Real> & normal) {

  UInt dim = normal.size();
  for (UInt i = 0; i < dim; ++i) {
    n[i] = normal[i] * tn;
    for (UInt j = 0; j < shapes.size(); ++j) {
      n[(1 + j) * dim + i] = -normal[i] * shapes[j];
    } 
  }
}

/* -------------------------------------------------------------------------- */
void Resolution::computeTalpha(Array<Real> & t_alpha, Vector<Real> & shapes,
			       Array<Real> & tangents) {
 
   for (auto && values:
	  zip(make_view(tangents, spatial_dimension),
	      make_view(t_alpha,  t_alpha.size()))) {

     auto & tangent = std::get<0>(values);
     auto & t       = std::get<1>(values);
     for (UInt i : arange(spatial_dimension)) {
       t[i] = -tangent[i];
       for (UInt j : arange(shapes.size())) {
	 t[(1 + j)*spatial_dimension + i] = -shapes[j] * tangent[i];
       }
    }
   }
}

/* -------------------------------------------------------------------------- */
void Resolution::computeNalpha(Array<Real> & n_alpha, Array<Real> & shapes_derivatives,
			       Vector<Real> & normal) {

  for (auto && values:
	 zip(make_view(shapes_derivatives, shapes_derivatives.size(),
		       n_alpha, n_alpha.size()))) {
    
    auto & shape_derivative = std::get<0>(values);
    auto & n                = std::get<1>(values);
    for (UInt i : arange(spatial_dimension)) {
      n[i] = 0;
      for (UInt j : arange(shapes.size())) {
	n[(1 + j)*spatial_dimension + i] = -shape_derivative[j]*normal[i];
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void Resolution::computeDalpha(Array<Real> & d_alpha, Array<Real> & n_alpha,
			       Array<Real> & t_alpha, Matrix<Real> & surface_matrix, Real gap) {


}

  
/* -------------------------------------------------------------------------- */
void Resolution::computeCoordiantes(const Element & el, Array<Real> & coords) {
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);
  Vector<UInt> connect = model.getMesh().getConnectivity(el.type, _not_ghost)
                                           .begin(nb_nodes_per_element)[el.element]; 

  for (UInt n = 0; n < nb_nodes_per_element; ++n) {
    UInt node = connect[n];
    for (UInt s: arange(spatial_dimension)) {
      coords(s, n) = this->positions(node, s);
    }
  }
}

