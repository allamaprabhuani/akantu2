/**
 * @file   solid_mechanics_model_boundary.cc
 * @author Guillaume ANCIAUX <guillaume.anciaux@epfl.ch>
 * @date   Fri Nov 19 10:23:03 2010
 *
 * @brief
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
#include "solid_mechanics_model.hh"
#include "dof_synchronizer.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
/**
 * @param myf pointer  to a function that fills a  vector/tensor with respect to
 * passed coordinates
 * @param function_type  flag to  specify the take  of function:  _bft_stress is
 * tensor like and _bft_forces is traction like.
 */
void SolidMechanicsModel::computeForcesFromFunction(BoundaryFunction myf,
						    BoundaryFunctionType function_type){
  /** function type is
   ** _bft_forces : traction function is given
   ** _bft_stress : stress function is given
   */
  GhostType ghost_type = _not_ghost;

  std::stringstream name; name << id << ":imposed_traction";
  Vector<Real> traction_funct(0, spatial_dimension, name.str());

  name.clear(); name << id << ":imposed_stresses";
  Vector<Real> stress_funct(0, spatial_dimension*spatial_dimension, name.str());

  name.clear(); name << id << ":quad_coords";
  Vector<Real> quad_coords(0, spatial_dimension, name.str());

  UInt offset = 0;
  switch(function_type) {
  case _bft_stress:
    offset = spatial_dimension * spatial_dimension; break;
  case _bft_forces:
    offset = spatial_dimension; break;
  }

  //prepare the loop over element types
  const Mesh::ConnectivityTypeList & type_list = getFEMBoundary().getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != getFEMBoundary().getElementDimension()) continue;

    UInt nb_quad    = getFEMBoundary().getNbQuadraturePoints(*it, ghost_type);
    UInt nb_element = getFEMBoundary().getMesh().getNbElement(*it, ghost_type);

    getFEMBoundary().interpolateOnQuadraturePoints(getFEMBoundary().getMesh().getNodes(),
						   quad_coords, spatial_dimension, *it, ghost_type);


    Real * imposed_val = NULL;
    switch(function_type) {
    case _bft_stress:
      stress_funct.resize(nb_element*nb_quad);
      imposed_val = stress_funct.values;
      break;
    case _bft_forces:
      traction_funct.resize(nb_element*nb_quad);
      imposed_val = traction_funct.values;
      break;
    }

    /// sigma/tractions on each quadrature points
    Real * qcoord = quad_coords.values;
    for (UInt el = 0; el < nb_element; ++el) {
      for (UInt q = 0; q < nb_quad; ++q) {
	myf(qcoord, imposed_val);
	imposed_val += offset;
	qcoord += spatial_dimension;
      }
    }

    switch(function_type) {
    case _bft_stress:
      computeForcesByStressTensor(stress_funct, *it, ghost_type); break;
    case _bft_forces:
      computeForcesByTractionVector(traction_funct, *it, ghost_type); break;
    }
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::computeForcesByStressTensor(const Vector<Real> & stresses,
						      const ElementType & type,
						      const GhostType & ghost_type){
  AKANTU_DEBUG_IN();

  UInt nb_element = getFEMBoundary().getMesh().getNbElement(type, ghost_type);
  UInt nb_quad = getFEMBoundary().getNbQuadraturePoints(type, ghost_type);

  // check dimension match
  AKANTU_DEBUG_ASSERT(Mesh::getSpatialDimension(type) == getFEMBoundary().getElementDimension(),
		      "element type dimension does not match the dimension of boundaries : " <<
		      getFEMBoundary().getElementDimension() << " != " <<
		      Mesh::getSpatialDimension(type));

  // check size of the vector
  AKANTU_DEBUG_ASSERT(stresses.getSize() == nb_quad*nb_element,
		      "the size of the vector should be the total number of quadrature points");

  // check number of components
  AKANTU_DEBUG_ASSERT(stresses.getNbComponent() == spatial_dimension*spatial_dimension,
		      "the number of components should be the dimension of 2-tensors");



  std::stringstream name;
  name << id << ":traction_boundary:" << type;
  Vector<Real> funct(nb_element*nb_quad, spatial_dimension,name.str());

  const Vector<Real> & normals_on_quad = getFEMBoundary().getNormalsOnQuadPoints(type, ghost_type);

  Math::matrix_vector(spatial_dimension, spatial_dimension,
		      stresses,normals_on_quad, funct);

  computeForcesByTractionVector(funct, type, ghost_type);
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::
computeForcesByTractionVector(const Vector<Real> & tractions,
			      const ElementType & type,
			      const GhostType & ghost_type){
  AKANTU_DEBUG_IN();

  UInt nb_element = getFEMBoundary().getMesh().getNbElement(type, ghost_type);
  UInt nb_nodes_per_element = getFEMBoundary().getMesh().getNbNodesPerElement(type);
  UInt nb_quad = getFEMBoundary().getNbQuadraturePoints(type, ghost_type);

  // check dimension match
  AKANTU_DEBUG_ASSERT(Mesh::getSpatialDimension(type)
		      == getFEMBoundary().getElementDimension(),
		      "element type dimension does not match "
		      <<"the dimension of boundaries : " <<
		      getFEMBoundary().getElementDimension() << " != " <<
		      Mesh::getSpatialDimension(type));

  // check size of the vector
  AKANTU_DEBUG_ASSERT(tractions.getSize() == nb_quad*nb_element,
		      "the size of the vector should be the "
		      << "total number of quadrature points");

  // check number of components
  AKANTU_DEBUG_ASSERT(tractions.getNbComponent() == spatial_dimension,
		      "the number of components should be "
		      << "the spatial dimension of the problem");


  // do a complete copy of the vector
  Vector<Real> funct(tractions, true);

  // extend the vector to multiply by the shapes (prepare to assembly)
  funct.extendComponentsInterlaced(nb_nodes_per_element,spatial_dimension);

  // multiply by the shapes
  Real * funct_val = funct.values;
  Real * shapes_val = (getFEMBoundary().getShapes(type)).values;
  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt q = 0; q < nb_quad; ++q) {
      for (UInt n = 0; n < nb_nodes_per_element; ++n,++shapes_val) {
	for (UInt i = 0; i < spatial_dimension; ++i) {
	  *funct_val++ *= *shapes_val;
	}
      }
    }
  }

  // allocate the vector that will contain the integrated values
  std::stringstream name;
  name << id << ":integral_boundary:" << type;
  Vector<Real> int_funct(nb_element, spatial_dimension*nb_nodes_per_element,name.str());

  //do the integration
  getFEMBoundary().integrate(funct, int_funct, spatial_dimension*nb_nodes_per_element, type, ghost_type);

  // assemble the result into force vector
  getFEMBoundary().assembleVector(int_funct,
				  const_cast<Vector<Real> &>(getForce()),
				  dof_synchronizer->getLocalDOFEquationNumbers(),
				  spatial_dimension,
				  type, ghost_type);
  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
