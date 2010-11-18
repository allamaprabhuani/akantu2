/**
 * @file   fem.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jul 16 11:03:02 2010
 *
 * @brief  Implementation of the FEM class
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
#include "fem.hh"
#include "mesh.hh"
#include "element_class.hh"
#include "aka_math.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
FEM::FEM(Mesh & mesh, UInt element_dimension, FEMID id, MemoryID memory_id) :
  Memory(memory_id), id(id) {
  AKANTU_DEBUG_IN();
  this->element_dimension = (element_dimension != 0) ?
    element_dimension : mesh.getSpatialDimension();

  init();

  this->mesh = &mesh;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FEM::init() {
  for(UInt t = _not_defined; t < _max_element_type; ++t) {
    this->shapes            [t] = NULL;
    this->shapes_derivatives[t] = NULL;
    this->jacobians         [t] = NULL;

    this->ghost_shapes            [t] = NULL;
    this->ghost_shapes_derivatives[t] = NULL;
    this->ghost_jacobians         [t] = NULL;
  }
}

/* -------------------------------------------------------------------------- */
FEM::~FEM() {
  AKANTU_DEBUG_IN();

  mesh = NULL;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FEM::initShapeFunctions(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  Real * coord = mesh->getNodes().values;
  UInt spatial_dimension = mesh->getSpatialDimension();

  const Mesh::ConnectivityTypeList & type_list = mesh->getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin();
      it != type_list.end();
      ++it) {

    ElementType type = *it;

    UInt element_type_spatial_dimension = Mesh::getSpatialDimension(type);
    UInt nb_nodes_per_element           = Mesh::getNbNodesPerElement(type);
    UInt size_of_shapes    = FEM::getShapeSize(type);
    UInt size_of_shapesd   = FEM::getShapeDerivativesSize(type);
    UInt nb_quadrature_points = FEM::getNbQuadraturePoints(type);

    if(element_type_spatial_dimension != element_dimension) continue;

    UInt * elem_val;
    UInt nb_element;
    std::string ghost = "";

    if(ghost_type == _not_ghost) {
      elem_val   = mesh->getConnectivity(type).values;
      nb_element = mesh->getConnectivity(type).getSize();
    } else {
      ghost = "ghost_";
      elem_val   = mesh->getGhostConnectivity(type).values;
      nb_element = mesh->getGhostConnectivity(type).getSize();
    }

    std::stringstream sstr_shapes;
    sstr_shapes << id << ":" << ghost << "shapes:" << type;
    Vector<Real> * shapes_tmp = &(alloc<Real>(sstr_shapes.str(),
					      nb_element*nb_quadrature_points,
					      size_of_shapes));

    std::stringstream sstr_shapesd;
    sstr_shapesd << id << ":" << ghost << "shapes_derivatives:" << type;
    Vector<Real> * shapes_derivatives_tmp = &(alloc<Real>(sstr_shapesd.str(),
							  nb_element*nb_quadrature_points,
							  size_of_shapesd));

    std::stringstream sstr_jacobians;
    sstr_jacobians << id << ":" << ghost << "jacobians:" << type;
    Vector<Real> * jacobians_tmp = &(alloc<Real>(sstr_jacobians.str(),
						 nb_element*nb_quadrature_points,
						 1));

    Real * shapes_val    = shapes_tmp->values;
    Real * shapesd_val   = shapes_derivatives_tmp->values;
    Real * jacobians_val = jacobians_tmp->values;

    /* -------------------------------------------------------------------------- */
    /* compute shapes when no rotation is required */

#define COMPUTE_SHAPES(type)						\
    do {								\
      Real local_coord[element_dimension * nb_nodes_per_element];	\
      for (UInt elem = 0; elem < nb_element; ++elem) {			\
	int offset = elem * nb_nodes_per_element;			\
	for (UInt id = 0; id < nb_nodes_per_element; ++id) {		\
	  memcpy(local_coord + id * spatial_dimension,			\
		 coord + elem_val[offset + id] * spatial_dimension,	\
		 spatial_dimension*sizeof(Real));			\
	}								\
	ElementClass<type>::preComputeStandards(local_coord,		\
						spatial_dimension,	\
						shapes_val,		\
						shapesd_val,		\
						jacobians_val);		\
	shapes_val += size_of_shapes*nb_quadrature_points;		\
	shapesd_val += size_of_shapesd*nb_quadrature_points;		\
	jacobians_val += nb_quadrature_points;				\
      }									\
    } while(0)

/* -------------------------------------------------------------------------- */

    switch(type) {
    case _line_1       : { COMPUTE_SHAPES(_line_1      ); break; }
    case _line_2       : { COMPUTE_SHAPES(_line_2      ); break; }
    case _triangle_1   : { COMPUTE_SHAPES(_triangle_1  ); break; }
    case _triangle_2   : { COMPUTE_SHAPES(_triangle_2  ); break; }
    case _tetrahedra_1 : { COMPUTE_SHAPES(_tetrahedra_1); break; }
    case _tetrahedra_2 : { COMPUTE_SHAPES(_tetrahedra_2); break; }
    case _point:
    case _not_defined:
    case _max_element_type:  {
      AKANTU_DEBUG_ERROR("Wrong type : " << type);
      break; }
    }
#undef COMPUTE_SHAPES

    if(ghost_type == _not_ghost) {
      shapes[type]             = shapes_tmp;
      shapes_derivatives[type] = shapes_derivatives_tmp;
      jacobians[type]          = jacobians_tmp;
    } else {
      ghost_shapes[type]             = shapes_tmp;
      ghost_shapes_derivatives[type] = shapes_derivatives_tmp;
      ghost_jacobians[type]          = jacobians_tmp;
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FEM::computeNormalsOnQuadPoints(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  Real * coord = mesh->getNodes().values;
  UInt spatial_dimension = mesh->getSpatialDimension();

  const Mesh::ConnectivityTypeList & type_list = mesh->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin();
      it != type_list.end();
      ++it) {

    ElementType type = *it;

    UInt element_type_spatial_dimension = Mesh::getSpatialDimension(type);
    UInt nb_nodes_per_element           = Mesh::getNbNodesPerElement(type);
    UInt nb_quad_points                 = FEM::getNbQuadraturePoints(type);

    if(element_type_spatial_dimension != element_dimension) continue;

    UInt * elem_val;
    UInt nb_element;
    std::string ghost = "";

    if(ghost_type == _not_ghost) {
      elem_val   = mesh->getConnectivity(type).values;
      nb_element = mesh->getConnectivity(type).getSize();
    } else {
      ghost = "ghost_";
      elem_val   = mesh->getGhostConnectivity(type).values;
      nb_element = mesh->getGhostConnectivity(type).getSize();
    }

    std::stringstream sstr_normals_on_quad;
    sstr_normals_on_quad << id << ":" << ghost << "normals_onquad:" << type;
    Vector<Real> * normals_on_quad_tmp = &(alloc<Real>(sstr_normals_on_quad.str(),
					      nb_element*nb_quad_points,
					      spatial_dimension));

    Real * normals_on_quad_val    = normals_on_quad_tmp->values;

#define COMPUTE_NORMALS_ON_QUAD(type)					\
    do {								\
      Real local_coord[spatial_dimension * nb_nodes_per_element];	\
      for (UInt elem = 0; elem < nb_element; ++elem) {			\
	int offset = elem * nb_nodes_per_element;			\
	for (UInt id = 0; id < nb_nodes_per_element; ++id) {		\
	  memcpy(local_coord + id * spatial_dimension,			\
		 coord + elem_val[offset + id] * spatial_dimension,	\
		 spatial_dimension*sizeof(Real));			\
	}								\
	ElementClass<type>::computeNormalsOnQuadPoint(local_coord,      \
						  spatial_dimension,	\
						  normals_on_quad_val);	\
	normals_on_quad_val += spatial_dimension*nb_quad_points;      	\
      }									\
    } while(0)

    switch(type) {
    case _line_1       : { COMPUTE_NORMALS_ON_QUAD(_line_1      ); break; }
    case _line_2       : { COMPUTE_NORMALS_ON_QUAD(_line_2      ); break; }
    case _triangle_1   : { COMPUTE_NORMALS_ON_QUAD(_triangle_1  ); break; }
    case _triangle_2   : { COMPUTE_NORMALS_ON_QUAD(_triangle_2  ); break; }
    case _tetrahedra_1 : { COMPUTE_NORMALS_ON_QUAD(_tetrahedra_1); break; }
    case _tetrahedra_2 : { COMPUTE_NORMALS_ON_QUAD(_tetrahedra_2); break; }
    case _point:
    case _not_defined:
    case _max_element_type:  {
      AKANTU_DEBUG_ERROR("Wrong type : " << type);
      break; }
    }
#undef COMPUTE_SHAPES

    if(ghost_type == _not_ghost) {
      normals_on_quad_points[type]             = normals_on_quad_tmp;
    } else {
      AKANTU_DEBUG_ERROR("to be implemented");
    }
  }
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void FEM::interpolateOnQuadraturePoints(const Vector<Real> &in_u,
					Vector<Real> &out_uq,
					UInt nb_degre_of_freedom,
					const ElementType & type,
					GhostType ghost_type,
					const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();

  Vector<Real> * shapes_loc;
  UInt nb_element;
  UInt * conn_val;

  if(ghost_type == _not_ghost) {
    shapes_loc = shapes[type];
    nb_element = mesh->getNbElement(type);
    conn_val   = mesh->getConnectivity(type).values;
  } else {
    shapes_loc = ghost_shapes[type];
    nb_element = mesh->getNbGhostElement(type);
    conn_val   = mesh->getGhostConnectivity(type).values;
  }

  AKANTU_DEBUG_ASSERT(shapes_loc != NULL,
		      "No shapes for the type " << type);

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = FEM::getNbQuadraturePoints(type);
  UInt size_of_shapes       = FEM::getShapeSize(type);

  AKANTU_DEBUG_ASSERT(in_u.getSize() == mesh->getNbNodes(),
		      "The vector in_u(" << in_u.getID()
		      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(in_u.getNbComponent() == nb_degre_of_freedom,
		      "The vector in_u(" << in_u.getID()
		      << ") has not the good number of component.");

  AKANTU_DEBUG_ASSERT(out_uq.getNbComponent() == nb_degre_of_freedom ,
		      "The vector out_uq(" << out_uq.getID()
		      << ") has not the good number of component.");

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element = filter_elements->getSize();
    filter_elem_val = filter_elements->values;
  }

  out_uq.resize(nb_element * nb_quadrature_points);

  Real * shape_val = shapes_loc->values;
  Real * u_val     = in_u.values;
  Real * uq_val    = out_uq.values;

  UInt offset_uq   = out_uq.getNbComponent()*nb_quadrature_points;

  Real * shape = shape_val;
  Real * u = static_cast<Real *>(calloc(nb_nodes_per_element * nb_degre_of_freedom,
					sizeof(Real)));

  for (UInt el = 0; el < nb_element; ++el) {
    UInt el_offset = el * nb_nodes_per_element;
    if(filter_elements != NULL) {
      shape     = shape_val + filter_elem_val[el] * size_of_shapes*nb_quadrature_points;
      el_offset = filter_elem_val[el] * nb_nodes_per_element;
    }

    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      memcpy(u + n * nb_degre_of_freedom,
	     u_val + conn_val[el_offset + n] * nb_degre_of_freedom,
	     nb_degre_of_freedom * sizeof(Real));
    }

    /// Uq = Shape * U : matrix product
    Math::matrix_matrix(nb_quadrature_points, nb_degre_of_freedom, nb_nodes_per_element,
			shape, u, uq_val);

    uq_val += offset_uq;
    if(filter_elements == NULL) {
      shape += size_of_shapes*nb_quadrature_points;
    }
  }

  free(u);

#undef INIT_VARIABLES
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FEM::gradientOnQuadraturePoints(const Vector<Real> &in_u,
				     Vector<Real> &out_nablauq,
				     UInt nb_degre_of_freedom,
				     const ElementType & type,
				     GhostType ghost_type,
				     const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();

  Vector<Real> * shapesd_loc;
  UInt nb_element;
  UInt * conn_val;

  if(ghost_type == _not_ghost) {
    shapesd_loc = shapes_derivatives[type];
    nb_element  = mesh->getNbElement(type);
    conn_val    = mesh->getConnectivity(type).values;
  } else {
    shapesd_loc = ghost_shapes_derivatives[type];
    nb_element  = mesh->getNbGhostElement(type);
    conn_val    = mesh->getGhostConnectivity(type).values;
  }

  AKANTU_DEBUG_ASSERT(shapesd_loc != NULL,
		      "No shapes for the type " << type);

  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(type);
  UInt size_of_shapes_derivatives = FEM::getShapeDerivativesSize(type);
  UInt nb_quadrature_points       = FEM::getNbQuadraturePoints(type);

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element = filter_elements->getSize();
    filter_elem_val = filter_elements->values;
  }

  AKANTU_DEBUG_ASSERT(in_u.getSize() == mesh->getNbNodes(),
		      "The vector in_u(" << in_u.getID()
		      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(in_u.getNbComponent() == nb_degre_of_freedom ,
		      "The vector in_u(" << in_u.getID()
		      << ") has not the good number of component.");

  AKANTU_DEBUG_ASSERT(out_nablauq.getNbComponent()
		      == nb_degre_of_freedom * element_dimension,
		      "The vector out_nablauq(" << out_nablauq.getID()
		      << ") has not the good number of component.");

  out_nablauq.resize(nb_element * nb_quadrature_points);

  Real * shaped_val  = shapesd_loc->values;
  Real * u_val       = in_u.values;
  Real * nablauq_val = out_nablauq.values;

  UInt offset_nablauq = nb_degre_of_freedom * element_dimension;
  UInt offset_shaped  = nb_nodes_per_element * element_dimension;

  Real * shaped  = shaped_val;
  Real * u       = static_cast<Real *>(calloc(nb_nodes_per_element * nb_degre_of_freedom,
					      sizeof(Real)));

  for (UInt el = 0; el < nb_element; ++el) {
    UInt el_offset = el * nb_nodes_per_element;
    if(filter_elements != NULL) {
      shaped    = shaped_val  + filter_elem_val[el] * size_of_shapes_derivatives*nb_quadrature_points;
      el_offset = filter_elem_val[el] * nb_nodes_per_element;
    }

    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      memcpy(u + n * nb_degre_of_freedom,
	     u_val + conn_val[el_offset + n] * nb_degre_of_freedom,
	     nb_degre_of_freedom * sizeof(Real));
    }

    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      /// \nabla(U) = U^t * dphi/dx
      Math::matrixt_matrix(nb_degre_of_freedom, element_dimension, nb_nodes_per_element,
			   u,
			   shaped,
			   nablauq_val);

      nablauq_val += offset_nablauq;
      shaped      += offset_shaped;
    }
  }

  free(u);

#undef INIT_VARIABLES
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FEM::integrate(const Vector<Real> & in_f,
		    Vector<Real> &intf,
		    UInt nb_degre_of_freedom,
		    const ElementType & type,
		    GhostType ghost_type,
		    const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();

  Vector<Real> * jac_loc;
  UInt nb_element;

  if(ghost_type == _not_ghost) {
    jac_loc     = jacobians[type];
    nb_element  = mesh->getNbElement(type);
  } else {
    jac_loc     = ghost_jacobians[type];
    nb_element  = mesh->getNbGhostElement(type);
  }

  UInt nb_quadrature_points = FEM::getNbQuadraturePoints(type);
  
  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element      = filter_elements->getSize();
    filter_elem_val = filter_elements->values;
  }

  AKANTU_DEBUG_ASSERT(in_f.getSize() == nb_element * nb_quadrature_points,
		      "The vector in_f(" << in_f.getID() << " size " << in_f.getSize()
		      << ") has not the good size (" << nb_element << ").");
  AKANTU_DEBUG_ASSERT(in_f.getNbComponent() == nb_degre_of_freedom ,
		      "The vector in_f(" << in_f.getID()
		      << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(intf.getNbComponent() == nb_degre_of_freedom,
		      "The vector intf(" << intf.getID()
		      << ") has not the good number of component.");


  intf.resize(nb_element);

  Real * in_f_val = in_f.values;
  Real * intf_val = intf.values;
  Real * jac_val  = jac_loc->values;

  UInt offset_in_f = in_f.getNbComponent()*nb_quadrature_points;
  UInt offset_intf = intf.getNbComponent();

  Real * jac      = jac_val;

  for (UInt el = 0; el < nb_element; ++el) {
    if(filter_elements != NULL) {
      jac      = jac_val  + filter_elem_val[el] * nb_quadrature_points;
    }

    integrate(in_f_val, jac, intf_val, nb_degre_of_freedom, nb_quadrature_points);

    in_f_val += offset_in_f;
    intf_val += offset_intf;
    if(filter_elements == NULL) {
      jac      += nb_quadrature_points;
    }
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
Real FEM::integrate(const Vector<Real> & in_f,
		    const ElementType & type,
		    GhostType ghost_type,
		    const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();
  Vector<Real> * jac_loc;
  UInt nb_element;

  if(ghost_type == _not_ghost) {
    jac_loc     = jacobians[type];
    nb_element  = mesh->getNbElement(type);
  } else {
    jac_loc     = ghost_jacobians[type];
    nb_element  = mesh->getNbGhostElement(type);
  }

  UInt nb_quadrature_points = FEM::getNbQuadraturePoints(type);
  
  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element      = filter_elements->getSize();
    filter_elem_val = filter_elements->values;
  }

  AKANTU_DEBUG_ASSERT(in_f.getSize() == nb_element * nb_quadrature_points,
		      "The vector in_f(" << in_f.getID()
		      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(in_f.getNbComponent() == 1,
		      "The vector in_f(" << in_f.getID()
		      << ") has not the good number of component.");

  Real intf = 0.;
  Real * in_f_val  = in_f.values;
  Real * jac_val   = jac_loc->values;
  UInt offset_in_f = in_f.getNbComponent();
  Real * jac       = jac_val;

  for (UInt el = 0; el < nb_element; ++el) {
    if(filter_elements != NULL) {
      jac = jac_val  + filter_elem_val[el] * nb_quadrature_points;
    }
    Real el_intf;
    integrate(in_f_val, jac, &el_intf, 1, nb_quadrature_points);
    intf += el_intf;

    in_f_val += offset_in_f;
    if(filter_elements == NULL) {
      jac += nb_quadrature_points;
    }
  }

  AKANTU_DEBUG_OUT();
  return intf;
}

/* -------------------------------------------------------------------------- */
void FEM::assembleVector(const Vector<Real> & elementary_vect,
			 Vector<Real> & nodal_values,
			 UInt nb_degre_of_freedom,
			 const ElementType & type,
			 GhostType ghost_type,
			 const Vector<UInt> * filter_elements,
			 Real scale_factor) const {
  AKANTU_DEBUG_IN();

  UInt nb_element;
  UInt * conn_val;

  if(ghost_type == _not_ghost) {
    nb_element  = mesh->getNbElement(type);
    conn_val    = mesh->getConnectivity(type).values;
  } else {
    nb_element  = mesh->getNbGhostElement(type);
    conn_val    = mesh->getGhostConnectivity(type).values;
  }

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_nodes = mesh->getNbNodes();

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element      = filter_elements->getSize();
    filter_elem_val = filter_elements->values;
  }

  AKANTU_DEBUG_ASSERT(elementary_vect.getSize() == nb_element,
		      "The vector elementary_vect(" << elementary_vect.getID()
		      << ") has not the good size.");

  AKANTU_DEBUG_ASSERT(elementary_vect.getNbComponent()
		      == nb_degre_of_freedom*nb_nodes_per_element,
		      "The vector elementary_vect(" << elementary_vect.getID()
		      << ") has not the good number of component.");

  AKANTU_DEBUG_ASSERT(nodal_values.getNbComponent() == nb_degre_of_freedom,
		      "The vector nodal_values(" << nodal_values.getID()
		      << ") has not the good number of component.");

  nodal_values.resize(nb_nodes);

  Real * elementary_vect_val = elementary_vect.values;
  Real * nodal_values_val    = nodal_values.values;

  for (UInt el = 0; el < nb_element; ++el) {
    UInt el_offset = el * nb_nodes_per_element;
    if(filter_elements != NULL) {
      el_offset = filter_elem_val[el] * nb_nodes_per_element;
    }
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt node = conn_val[el_offset + n];
      UInt offset_node = node * nb_degre_of_freedom;
      for (UInt d = 0; d < nb_degre_of_freedom; ++d) {
	nodal_values_val[offset_node + d] += scale_factor * elementary_vect_val[d];
      }
      elementary_vect_val += nb_degre_of_freedom;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FEM::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "FEM [" << std::endl;
  stream << space << " + id                : " << id << std::endl;
  stream << space << " + element dimension : " << element_dimension << std::endl;

  stream << space << " + mesh [" << std::endl;
  mesh->printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << " + connectivity type information [" << std::endl;
  const Mesh::ConnectivityTypeList & type_list = mesh->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if (mesh->getSpatialDimension(*it) != element_dimension) continue;
    stream << space << AKANTU_INDENT << AKANTU_INDENT << " + " << *it <<" [" << std::endl;
    if(shapes[*it]) {
      shapes            [*it]->printself(stream, indent + 3);
      shapes_derivatives[*it]->printself(stream, indent + 3);
      jacobians         [*it]->printself(stream, indent + 3);
    }
    stream << space << AKANTU_INDENT << AKANTU_INDENT << "]" << std::endl;
  }
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << "]" << std::endl;
}


__END_AKANTU__
