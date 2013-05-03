/**
 * @file   solid_mechanics_model_boundary.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Nov 19 15:27:01 2010
 *
 * @brief  implementation of functions to impose boundary conditions
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
#include "aka_types.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

// XXX TODO FIXME delete this
//class WarppingSurfaceLoadFunctor : public SolidMechanicsModel::SurfaceLoadFunctor {
//public:
//  WarppingSurfaceLoadFunctor(BoundaryFunction function) : function(function) {}
//
//  void traction(const Vector<Real> & position,
//		Vector<Real> & force,
//		const Vector<Real> & normal,
//		Surface surface_id) {
//    function(position.storage(),
//	     force.storage(),
//	     normal.storage(),
//	     surface_id);
//  }
//  void stress(const Vector<Real> & position,
//	      Matrix<Real> & stress,
//	      const Vector<Real> & normal,
//	      Surface surface_id) {
//    function(position.storage(),
//	     stress.storage(),
//	     normal.storage(),
//	     surface_id);
//  }
//private:
//  BoundaryFunction function;
//};

/**
 * @param myf pointer  to a function that fills a  vector/tensor with respect to
 * passed coordinates
 * @param function_type  flag to  specify the take  of function:  _bft_stress is
 * tensor like and _bft_forces is traction like.
 */
//void SolidMechanicsModel::computeForcesFromFunction(BoundaryFunction function,
//						    BoundaryFunctionType function_type){
//  WarppingSurfaceLoadFunctor functor(function);
//  computeForcesFromFunction(functor, function_type);
//}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModel::computeForcesByStressTensor(const Array<Real> & stresses,
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
  Array<Real> funct(nb_element*nb_quad, spatial_dimension,name.str());

  const Array<Real> & normals_on_quad = getFEMBoundary().getNormalsOnQuadPoints(type, ghost_type);

  Math::matrix_vector(spatial_dimension, spatial_dimension,
		      stresses,normals_on_quad, funct);

  computeForcesByTractionArray(funct, type, ghost_type);
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::
computeForcesByTractionArray(const Array<Real> & tractions,
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
  Array<Real> funct(tractions, true);

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
  Array<Real> int_funct(nb_element, spatial_dimension*nb_nodes_per_element,name.str());

  //do the integration
  getFEMBoundary().integrate(funct, int_funct, spatial_dimension*nb_nodes_per_element, type, ghost_type);

  // assemble the result into force vector
  getFEMBoundary().assembleArray(int_funct,
				  *force,
				  dof_synchronizer->getLocalDOFEquationNumbers(),
				  spatial_dimension,
				  type, ghost_type);
  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
