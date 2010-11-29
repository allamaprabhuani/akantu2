/**
 * @file   solid_mechanics_model_boundary.cc
 * @author Guillaume ANCIAUX <anciaux@lsmscluster1.epfl.ch>
 * @date   Fri Nov 19 10:23:03 2010
 *
 * @brief  
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
/** 
 * @param myf pointer to a function that fills a vector/tensor with respect to passed coordinates
 * @param function_type flag to specify the take of function: 0 is tensor like and 1 is traction like. 
 */
void SolidMechanicsModel::computeForcesFromFunction(void (*myf)(double *,double *), 
						    BoundaryFunctionType function_type){
  

  /** function type is 
   ** 0 : stress function is given
   ** 1 : traction function is given
   */

  std::stringstream name;
  name << id << ":solidmechanics:imposed_traction";
  Vector<Real> traction_funct(0, spatial_dimension,name.str());
  name.clear();
  name << id << ":solidmechanics:imposed_stresses";  
  Vector<Real> stress_funct(0, spatial_dimension*spatial_dimension,name.str());
  name.clear();
  name << id << ":solidmechanics:quad_coords";  
  Vector<Real> quad_coords(0,spatial_dimension,"quad_coords");

  UInt offset = 0;
  switch(function_type) {
  case _bft_stress:
    offset = spatial_dimension*spatial_dimension; break;
  case _bft_forces:
    offset = spatial_dimension; break;
  }


  //prepare the loop over element types
  const Mesh::ConnectivityTypeList & type_list = fem_boundary->getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != fem_boundary->getElementDimension()) continue;

    UInt nb_quad              = FEM::getNbQuadraturePoints(*it);
    UInt nb_element = fem_boundary->getMesh().getNbElement(*it);

    fem_boundary->interpolateOnQuadraturePoints(fem_boundary->getMesh().getNodes(), 
						quad_coords, spatial_dimension, (*it));


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
    for (UInt el = 0; el < nb_element; ++el) {
      Real * qcoord = quad_coords.values+el*nb_quad*spatial_dimension;
      for (UInt q = 0; q < nb_quad; ++q) {
	myf(qcoord+q,imposed_val);
	imposed_val += offset;
      }
    }

    switch(function_type) {
    case _bft_stress:
      computeForcesByStressTensor(stress_funct,(*it)); break;
    case _bft_forces:
      computeForcesByTractionVector(traction_funct,(*it)); break;
    }
  }
} 

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::computeForcesByStressTensor(const Vector<Real> & stresses, const ElementType & type){
  AKANTU_DEBUG_IN();

  UInt nb_element = fem_boundary->getMesh().getNbElement(type);
  UInt nb_quad = FEM::getNbQuadraturePoints(type);

  // check dimension match
  AKANTU_DEBUG_ASSERT(Mesh::getSpatialDimension(type) == fem_boundary->getElementDimension(),
		      "element type dimension does not match the dimension of boundaries : " << 
		      fem_boundary->getElementDimension() << " != " <<
		      Mesh::getSpatialDimension(type));

  // check size of the vector
  AKANTU_DEBUG_ASSERT(stresses.getSize() == nb_quad*nb_element,
		      "the size of the vector should be the total number of quadrature points");

  // check number of components
  AKANTU_DEBUG_ASSERT(stresses.getNbComponent() == spatial_dimension*spatial_dimension,
		      "the number of components should be the dimension of 2-tensors");



  std::stringstream name;
  name << id << ":solidmechanics:" << type << ":traction_boundary";
  Vector<Real> funct(nb_element*nb_quad, spatial_dimension,name.str());
  const Vector<Real> & normals_on_quad = fem_boundary->getNormalsOnQuadPoints(type);  

  Math::matrix_vector(spatial_dimension,spatial_dimension,stresses,normals_on_quad,funct);
  computeForcesByTractionVector(funct,type);
  AKANTU_DEBUG_OUT();  
}
/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::computeForcesByTractionVector(const Vector<Real> & tractions, const ElementType & type){
  AKANTU_DEBUG_IN();

  UInt nb_element = fem_boundary->getMesh().getNbElement(type);
  UInt nb_nodes_per_element = fem_boundary->getMesh().getNbNodesPerElement(type);
  UInt nb_quad = FEM::getNbQuadraturePoints(type);

  // check dimension match
  AKANTU_DEBUG_ASSERT(Mesh::getSpatialDimension(type) == fem_boundary->getElementDimension(),
		      "element type dimension does not match the dimension of boundaries : " << 
		      fem_boundary->getElementDimension() << " != " <<
		      Mesh::getSpatialDimension(type));

  // check size of the vector
  AKANTU_DEBUG_ASSERT(tractions.getSize() == nb_quad*nb_element,
		      "the size of the vector should be the total number of quadrature points");

  // check number of components
  AKANTU_DEBUG_ASSERT(tractions.getNbComponent() == spatial_dimension,
		      "the number of components should be the spatial dimension of the problem");


  // do a complete copy of the vector
  Vector<Real> funct(tractions,true);
  // extend the vector to multiply by the shapes (prepare to assembly)
  funct.extendComponentsInterlaced(nb_nodes_per_element,spatial_dimension);
  // multiply by the shapes
  Real * funct_val = funct.values;
  Real * shapes_val = (fem_boundary->getShapes(type)).values;
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
  name << id << ":solidmechanics:" << type << ":integral_boundary";
  Vector<Real> int_funct(nb_element, spatial_dimension*nb_nodes_per_element,name.str());
  //do the integration
  fem_boundary->integrate(funct, int_funct, spatial_dimension*nb_nodes_per_element, type);
  // assemble the result into force vector
  fem_boundary->assembleVector(int_funct,const_cast<Vector<Real> &>(getForce()), spatial_dimension, type);
  AKANTU_DEBUG_OUT();  
}

__END_AKANTU__
