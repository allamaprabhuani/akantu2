/**
 * @file   fem.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jul 16 11:03:02 2010
 *
 * @brief  Implementation of the FEM class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "fem.hh"
#include "mesh.hh"
#include "element_class.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
FEM::FEM(UInt spatial_dimension, FEMID id, MemoryID memory_id) :
  Memory(memory_id), id(id), spatial_dimension(spatial_dimension), created_mesh(true) {
  AKANTU_DEBUG_IN();
  std::stringstream sstr;
  sstr << id << ":mesh";

  for(UInt t = _not_defined; t < _max_element_type; ++t) {
    shapes[t] = NULL;
    shapes_derivatives[t] = NULL;
    jacobians[t] = NULL;
  }

  this->mesh = new Mesh(spatial_dimension, sstr.str(), memory_id);
  AKANTU_DEBUG_OUT();
};

/* -------------------------------------------------------------------------- */
FEM::FEM(Mesh & mesh, UInt spatial_dimension, FEMID id, MemoryID memory_id) :
  Memory(memory_id), id(id), created_mesh(false) {
  AKANTU_DEBUG_IN();
  this->spatial_dimension = (spatial_dimension != 0) ?
    spatial_dimension : mesh.getSpatialDimension();

  for(UInt t = _not_defined; t < _max_element_type; ++t) {
    shapes[t] = NULL;
    shapes_derivatives[t] = NULL;
    jacobians[t] = NULL;
  }

  this->mesh = &mesh;
  AKANTU_DEBUG_OUT();
};

/* -------------------------------------------------------------------------- */
FEM::~FEM() {
  AKANTU_DEBUG_IN();

  const Mesh::TypeList & type_list = mesh->getTypeList();
  Mesh::TypeList::const_iterator it;

  for(it = type_list.begin();
      it != type_list.end();
      ++it) {

    AKANTU_DEBUG(dblAccessory, "Deleting shapes vector of type " << *it);
    dealloc(shapes[*it]->getID());
    shapes[*it] = NULL;

    AKANTU_DEBUG(dblAccessory, "Deleting shapes derivatives vector of type " << *it);
    dealloc(shapes_derivatives[*it]->getID());
    shapes_derivatives[*it] = NULL;

    AKANTU_DEBUG(dblAccessory, "Deleting jacobians vector of type " << *it);
    dealloc(jacobians[*it]->getID());
    jacobians[*it] = NULL;
  }

  if(created_mesh) {
    AKANTU_DEBUG(dblAccessory, "Deleting mesh");
    delete mesh;
  }

  mesh = NULL;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FEM::initShapeFunctions() {
  AKANTU_DEBUG_IN();
  Real * coord = mesh->getNodes().values;

  const Mesh::TypeList & type_list = mesh->getTypeList();
  Mesh::TypeList::const_iterator it;

  for(it = type_list.begin();
      it != type_list.end();
      ++it) {

    ElementType type = *it;

    UInt nb_nodes_per_element;
    UInt size_of_shapes;
    UInt size_of_shapesd;
    UInt size_of_jacobians;

#define INIT_VARIABLES(type)						\
    do {								\
      if(ElementClass<type>::getSpatialDimension() != spatial_dimension) \
	continue;							\
      nb_nodes_per_element =						\
	ElementClass<type>::getNbNodesPerElement();			\
      size_of_shapes =							\
	ElementClass<type>::getShapeSize();				\
      size_of_shapesd =							\
	ElementClass<type>::getShapeDerivatiesSize();			\
      size_of_jacobians =						\
	ElementClass<type>::getJacobiansSize();				\
    } while(0)

    switch(type) {
    case _line_1       : { INIT_VARIABLES(_line_1      ); break; }
    case _line_2       : { INIT_VARIABLES(_line_2      ); break; }
    case _triangle_1   : { INIT_VARIABLES(_triangle_1  ); break; }
    case _triangle_2   : { INIT_VARIABLES(_triangle_2  ); break; }
    case _tetrahedra_1 : { INIT_VARIABLES(_tetrahedra_1); break; }
    case _tetrahedra_2 : { INIT_VARIABLES(_tetrahedra_2); break; }
    case _not_defined:
    case _max_element_type:  {
      AKANTU_DEBUG_ERROR("Wrong type : " << type);
      break; }
    }

    UInt * elem_val  = mesh->getConnectivity(type).values;
    UInt nb_element = mesh->getConnectivity(type).getSize();

    std::stringstream sstr_shapes;
    sstr_shapes << id << ":shapes:" << type;
    shapes[type] = &(alloc<Real>(sstr_shapes.str(),
				 nb_element,
				 size_of_shapes));

    std::stringstream sstr_shapesd;
    sstr_shapesd << id << ":shapes_derivatives:" << type;
    shapes_derivatives[type] = &(alloc<Real>(sstr_shapesd.str(),
					     nb_element,
					     size_of_shapesd));

    std::stringstream sstr_jacobians;
    sstr_jacobians << id << ":jacobians:" << type;
    jacobians[type] = &(alloc<Real>(sstr_jacobians.str(),
				    nb_element,
				    size_of_jacobians));

    Real * shapes_val    = shapes[type]->values;
    Real * shapesd_val   = shapes_derivatives[type]->values;
    Real * jacobians_val = jacobians[type]->values;

#define COMPUTE_SHAPES(type)						\
    do {								\
      Real local_coord[spatial_dimension * nb_nodes_per_element];	\
      for (UInt elem = 0; elem < nb_element; ++elem) {			\
	int offset = elem * nb_nodes_per_element;			\
	for (UInt id = 0; id < nb_nodes_per_element; ++id) {		\
	  memcpy(local_coord + id * spatial_dimension,			\
		 coord + elem_val[offset + id] * spatial_dimension,	\
		 spatial_dimension*sizeof(Real));			\
	}								\
	ElementClass<type>::shapeFunctions(local_coord,			\
					   shapes_val,			\
					   shapesd_val,			\
					   jacobians_val);		\
	shapes_val += size_of_shapes;					\
	shapesd_val += size_of_shapesd;					\
	jacobians_val += size_of_jacobians;				\
      }									\
    } while(0)

    switch(type) {
    case _line_1       : { COMPUTE_SHAPES(_line_1      ); break; }
    case _line_2       : { COMPUTE_SHAPES(_line_2      ); break; }
    case _triangle_1   : { COMPUTE_SHAPES(_triangle_1  ); break; }
    case _triangle_2   : { COMPUTE_SHAPES(_triangle_2  ); break; }
    case _tetrahedra_1 : { COMPUTE_SHAPES(_tetrahedra_1); break; }
    case _tetrahedra_2 : { COMPUTE_SHAPES(_tetrahedra_2); break; }
    case _not_defined:
    case _max_element_type:  {
      AKANTU_DEBUG_ERROR("Wrong type : " << type);
      break; }
    }
  }
  AKANTU_DEBUG_OUT();
#undef INIT_VARIABLES
#undef COMPUTE_SHAPES
};

/* -------------------------------------------------------------------------- */
void FEM::interpolateOnQuadraturePoints(const Vector<Real> &inval,
					Vector<Real> &valonquad,
					ElementType type,
					const Vector<UInt> * elements){
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(shapes[type] != NULL,
		      "No shapes for the type " << type);


  UInt nb_nodes_per_element;
  UInt size_of_shapes;
  UInt nb_quadrature_points;

#define INIT_VARIABLES(type)						\
  do {									\
    nb_nodes_per_element =						\
      ElementClass<type>::getNbNodesPerElement();			\
    size_of_shapes =							\
      ElementClass<type>::getShapeSize();				\
    nb_quadrature_points =						\
      ElementClass<type>::getNbQuadraturePoints();			\
  } while(0)

  switch(type) {
  case _line_1       : { INIT_VARIABLES(_line_1      ); break; }
  case _line_2       : { INIT_VARIABLES(_line_2      ); break; }
  case _triangle_1   : { INIT_VARIABLES(_triangle_1  ); break; }
  case _triangle_2   : { INIT_VARIABLES(_triangle_2  ); break; }
  case _tetrahedra_1 : { INIT_VARIABLES(_tetrahedra_1); break; }
  case _tetrahedra_2 : { INIT_VARIABLES(_tetrahedra_2); break; }
  case _not_defined:
  case _max_element_type:  {
    AKANTU_DEBUG_ERROR("Wrong type : " << type);
    break; }
  }

  UInt * conn_val = mesh->getConnectivity(type).values;
  UInt nb_element = mesh->getConnectivity(type).getSize();

  UInt * elem_val = NULL;
  if(elements != NULL) {
    nb_element = elements->getSize();
    elem_val = elements->values;
  }
  UInt degre_of_freedom = inval.getNbComponent();
  valonquad.resize(nb_element * nb_quadrature_points);

  Real * shape_val = shapes[type]->values;
  Real * u_val  = inval.values;
  Real * uq_val = valonquad.values;

  UInt offset_shape = size_of_shapes;
  UInt offset_uq    = degre_of_freedom * nb_quadrature_points;

  Real * shape = shape_val;
  Real * uq    = uq_val;
  Real * u = static_cast<Real *>(calloc(nb_nodes_per_element * degre_of_freedom,
					sizeof(Real)));

  for (UInt el = 0; el < nb_element; ++el) {
    UInt el_offset = el * nb_nodes_per_element;
    if(elements != NULL) {
      shape     = shape_val + elem_val[el] * offset_shape;
      uq        = uq_val    + elem_val[el] * offset_uq;
      el_offset = elem_val[el] * nb_nodes_per_element;
    }

    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      memcpy(u + n * degre_of_freedom,
	     u_val + conn_val[el_offset + n] * degre_of_freedom,
	     degre_of_freedom * sizeof(Real));
    }

    /// Uq = Shape * U : matrix product
#ifdef AKANTU_USE_BLAS
#else
    for (UInt i = 0; i < nb_quadrature_points; ++i) {
      UInt uq_i = i * degre_of_freedom;
      UInt sh_i = i * nb_nodes_per_element;
      for (UInt j = 0; j < degre_of_freedom; ++j) {
	uq[uq_i + j] = 0.;
	for (UInt k = 0; k < nb_nodes_per_element; ++k) {
	  uq[uq_i + j] += shape[sh_i + k] * u[k * degre_of_freedom + j];
	}
      }
    }
#endif
    if(elements == NULL) {
      shape += offset_shape;
      uq    += offset_uq;
    }
  }

  free(u);

#undef INIT_VARIABLES
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FEM::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "FEM [" << std::endl;
  stream << space << " + id                : " << id << std::endl;
  stream << space << " + spatial dimension : " << spatial_dimension << std::endl;
  stream << space << " + mesh [" << std::endl;
  mesh->printself(stream, indent + 2);
  stream << space << " ]" << std::endl;
  stream << space << "]" << std::endl;
};


__END_AKANTU__
