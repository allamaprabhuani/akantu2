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
#include "sparse_matrix.hh"
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
  
  auto & internal_force =
    const_cast<Array<Real> &>(model.getInternalForce());
  
  const auto local_nodes =
    model.getMesh().getElementGroup(name).getNodes();

  auto & nodal_area =
    const_cast<Array<Real> &>(model.getNodalArea());
    
  auto & contact_map = model.getContactMap();
  
  for (auto & slave: local_nodes) {

    if (contact_map.find(slave) == contact_map.end())
      continue; 
    
    auto & element = contact_map[slave];

    const auto & conn = element.connectivity;
    const auto & type = element.master.type;

    auto nb_nodes_master = Mesh::getNbNodesPerElement(type);

    Vector<Real> shapes(nb_nodes_master);
    Matrix<Real> dnds(spatial_dimension - 1, nb_nodes_master);
       
#define GET_SHAPES_NATURAL(type)				\
    ElementClass<type>::computeShapes(element.projection, shapes)
    AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPES_NATURAL);
#undef GET_SHAPES_NATURAL  

#define GET_SHAPE_DERIVATIVES_NATURAL(type)				\
    ElementClass<type>::computeDNDS(element.projection, dnds)
    AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_DERIVATIVES_NATURAL);
#undef GET_SHAPE_DERIVATIVES_NATURAL
        
    Vector<Real> fc(conn.size() * spatial_dimension);
   
    Matrix<Real> surface_matrix(spatial_dimension - 1, spatial_dimension - 1);
    computeSurfaceMatrix(element.tangents, surface_matrix);

    Vector<Real> n(conn.size() * spatial_dimension);
    
    computeN(n, shapes, element.normal);

    computeNormalForce(fc, n, element.gap);
    
    Array<Real> t_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);
    Array<Real> n_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);
    Array<Real> d_alpha(conn.size() * spatial_dimension, spatial_dimension - 1);
   
    computeTalpha(t_alpha, shapes,             element.tangents);
    computeNalpha(n_alpha, dnds,               element.normal);
    computeDalpha(d_alpha, n_alpha, t_alpha, surface_matrix, element.gap);

    //computeFrictionForce(fc, d_alpha, gap);

    UInt nb_degree_of_freedom = internal_force.getNbComponent();
    for (UInt i = 0; i < conn.size(); ++i) {

      UInt n = conn[i];
      for (UInt j = 0; j < nb_degree_of_freedom; ++j) {	
	UInt offset_node = n * nb_degree_of_freedom + j;
	internal_force[offset_node] += fc[i*nb_degree_of_freedom + j];
	internal_force[offset_node] *= nodal_area[n];
      }
    }    
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Resolution::assembleStiffnessMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & stiffness =
    const_cast<SparseMatrix &>(model.getDOFManager().getMatrix("K"));

  const auto local_nodes =
    model.getMesh().getElementGroup(name).getNodes();

  auto & nodal_area =
    const_cast<Array<Real> &>(model.getNodalArea());
  
  auto & contact_map = model.getContactMap();
  
  for (auto & slave: local_nodes) {

    if (contact_map.find(slave) == contact_map.end()) {
      continue;
    }

    auto & master       = contact_map[slave].master;
    auto & gap          = contact_map[slave].gap;
    auto & projection   = contact_map[slave].projection;
    auto & normal       = contact_map[slave].normal;
    const auto & connectivity = contact_map[slave].connectivity;
    const ElementType & type  = master.type;

    UInt nb_nodes_master = Mesh::getNbNodesPerElement(master.type);

    Vector<Real> shapes(nb_nodes_master);
    Matrix<Real> shapes_derivatives(spatial_dimension - 1, nb_nodes_master);
       
#define GET_SHAPES_NATURAL(type)				\
    ElementClass<type>::computeShapes(projection, shapes)
    AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPES_NATURAL);
#undef GET_SHAPES_NATURAL  

#define GET_SHAPE_DERIVATIVES_NATURAL(type)				\
    ElementClass<type>::computeDNDS(projection, shapes_derivatives)
    AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_DERIVATIVES_NATURAL);
#undef GET_SHAPE_DERIVATIVES_NATURAL

    Matrix<Real> elementary_stiffness(connectivity.size() * spatial_dimension,
				      connectivity.size() * spatial_dimension);

    Matrix<Real> tangents(spatial_dimension - 1, spatial_dimension);
    Matrix<Real> global_coords(nb_nodes_master, spatial_dimension);

    computeCoordinates(master, global_coords);
    computeTangents(shapes_derivatives, global_coords, tangents);
    
    Matrix<Real> surface_matrix(spatial_dimension - 1, spatial_dimension - 1);
    computeSurfaceMatrix(tangents, surface_matrix);
    
    Vector<Real> n(connectivity.size() * spatial_dimension);
    Array<Real> t_alpha(connectivity.size() * spatial_dimension, spatial_dimension - 1);
    Array<Real> n_alpha(connectivity.size() * spatial_dimension, spatial_dimension - 1);
    Array<Real> d_alpha(connectivity.size() * spatial_dimension, spatial_dimension - 1);

    computeN(      n,        shapes,             normal);
    computeTalpha( t_alpha,  shapes,             tangents);
    computeNalpha( n_alpha,  shapes_derivatives, normal);
    computeDalpha( d_alpha,  n_alpha,  t_alpha,  surface_matrix, gap);
    
    Matrix<Real> kc(connectivity.size() * spatial_dimension,
		    connectivity.size() * spatial_dimension);
    computeTangentModuli(kc, n, n_alpha, d_alpha, surface_matrix, gap);
        
    std::vector<UInt> equations;
    UInt nb_degree_of_freedom = model.getSpatialDimension();

    std::vector<Real> areas;
    for (UInt i : arange(connectivity.size())) {
      UInt n = connectivity[i];
      for (UInt j : arange(nb_degree_of_freedom)) {
	equations.push_back(n * nb_degree_of_freedom + j);
	areas.push_back(nodal_area[n]);
      }
    }

    for (UInt i : arange(kc.rows())) {
      UInt row = equations[i];
      for (UInt j : arange(kc.cols())) {
	UInt col = equations[j];
	kc(i, j) *= areas[i];
	stiffness.add(row, col, kc(i, j));
      }	
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Resolution::computeTangents(Matrix<Real> & shapes_derivatives, Matrix<Real> & global_coords,
				 Matrix<Real> & tangents) {

  tangents.mul<false, false>(shapes_derivatives, global_coords);
}

/* -------------------------------------------------------------------------- */
void Resolution::computeSurfaceMatrix(Matrix<Real> & tangents, Matrix<Real> & surface_matrix) {

  surface_matrix.mul<false, true>(tangents, tangents);
  surface_matrix = surface_matrix.inverse();
}

/* -------------------------------------------------------------------------- */
void Resolution::computeN(Vector<Real> & n, Vector<Real> & shapes, Vector<Real> & normal) {

  UInt dim = normal.size();
  for (UInt i = 0; i < dim; ++i) {
    n[i] = normal[i];
    for (UInt j = 0; j < shapes.size(); ++j) {
      n[(1 + j) * dim + i] = -normal[i] * shapes[j];
    } 
  }
}

/* -------------------------------------------------------------------------- */
void Resolution::computeTalpha(Array<Real> & t_alpha, Vector<Real> & shapes,
			       Matrix<Real> & tangents) {
  t_alpha.clear();
  for (auto && values:
	  zip(tangents.transpose(),
	      make_view(t_alpha, t_alpha.size()))) {

     auto & tangent = std::get<0>(values);
     auto & t_s     = std::get<1>(values);
     for (UInt i : arange(spatial_dimension)) {
       t_s[i] = -tangent(i);
       for (UInt j : arange(shapes.size())) {
	 t_s[(1 + j)*spatial_dimension + i] = -shapes[j] * tangent(i);
       }
     }
  }
}

/* -------------------------------------------------------------------------- */
void Resolution::computeNalpha(Array<Real> & n_alpha, Matrix<Real> & shapes_derivatives,
			       Vector<Real> & normal) {
  n_alpha.clear();
  for (auto && values:
	 zip(shapes_derivatives.transpose(),
	     make_view(n_alpha, n_alpha.size()))) {
    auto & dnds  = std::get<0>(values);       
    auto & n_s   = std::get<1>(values);
    for (UInt i : arange(spatial_dimension)) {
      n_s[i] = 0;
      for (UInt j : arange(dnds.size())) {
	n_s[(1 + j)*spatial_dimension + i] = -dnds(j)*normal[i];
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void Resolution::computeDalpha(Array<Real> & d_alpha, Array<Real> & n_alpha,
			       Array<Real> & t_alpha, Matrix<Real> & surface_matrix,
			       Real & gap) {

  d_alpha.clear();
  for (auto && entry : zip(surface_matrix.transpose(),
			   make_view(d_alpha, d_alpha.size()))) {
    auto & a_s = std::get<0>(entry);
    auto & d_s = std::get<1>(entry);
    for (auto && values :
	   zip(arange(t_alpha.size()),
	       make_view(t_alpha, t_alpha.size()),
	       make_view(n_alpha, n_alpha.size()))) {
      auto & index = std::get<0>(values);
      auto & t_s   = std::get<1>(values);
      auto & n_s   = std::get<2>(values);
      
      d_s += (t_s + gap  * n_s);
      d_s *= a_s(index);
    }
  }
}

  
/* -------------------------------------------------------------------------- */
void Resolution::computeCoordinates(const Element & el, Matrix<Real> & coords) {
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);
  Vector<UInt> connect = model.getMesh().getConnectivity(el.type, _not_ghost)
    .begin(nb_nodes_per_element)[el.element]; 

  // change this to current position
  auto & positions = model.getMesh().getNodes();
  for (UInt n = 0; n < nb_nodes_per_element; ++n) {
    UInt node = connect[n];
    for (UInt s: arange(spatial_dimension)) {
      coords(n, s) = positions(node, s);
    }
  }
}

} // akantu
