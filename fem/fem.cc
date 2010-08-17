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
#include "aka_math.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
FEM::FEM(UInt spatial_dimension, FEMID id, MemoryID memory_id) :
  Memory(memory_id), id(id), spatial_dimension(spatial_dimension), created_mesh(true) {
  AKANTU_DEBUG_IN();
  std::stringstream sstr;
  sstr << id << ":mesh";

  for(UInt t = _not_defined; t < _max_element_type; ++t) {
    this->shapes[t]             = NULL;
    this->shapes_derivatives[t] = NULL;
    this->jacobians[t]          = NULL;
  }

  this->mesh = new Mesh(spatial_dimension, sstr.str(), memory_id);
  AKANTU_DEBUG_OUT();
}

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
}

/* -------------------------------------------------------------------------- */
FEM::~FEM() {
  AKANTU_DEBUG_IN();

  const Mesh::ConnectivityTypeList & type_list = mesh->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

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

  const Mesh::ConnectivityTypeList & type_list = mesh->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin();
      it != type_list.end();
      ++it) {

    ElementType type = *it;

    UInt nb_nodes_per_element = 0;
    UInt element_type_spatial_dimension = 0;
    UInt size_of_shapes = 0;
    UInt size_of_shapesd = 0;
    UInt size_of_jacobians = 0;

#define INIT_VARIABLES(type)						\
    do {								\
      element_type_spatial_dimension =					\
	ElementClass<type>::getSpatialDimension();			\
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

    if(element_type_spatial_dimension != spatial_dimension) continue;

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
}

/* -------------------------------------------------------------------------- */
void FEM::interpolateOnQuadraturePoints(const Vector<Real> &in_u,
					Vector<Real> &out_uq,
					UInt nb_degre_of_freedom,
					const ElementType & type,
					const Vector<UInt> * filter_elements){
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(shapes[type] != NULL,
		      "No shapes for the type " << type);

  UInt nb_nodes_per_element = mesh->getConnectivity(type).getNbComponent();
  UInt size_of_shapes = shapes[type]->getNbComponent();
  UInt nb_quadrature_points = getNbQuadraturePoints(type);
  UInt nb_element = mesh->getConnectivity(type).getSize();

  AKANTU_DEBUG_ASSERT(in_u.getSize() == mesh->getNbNodes(),
		      "The vector in_u(" << in_u.getID()
		      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(in_u.getNbComponent() == nb_degre_of_freedom,
		      "The vector in_u(" << in_u.getID()
		      << ") has not the good number of component.");

  AKANTU_DEBUG_ASSERT(out_uq.getNbComponent() == nb_degre_of_freedom * nb_quadrature_points,
		      "The vector out_uq(" << out_uq.getID()
		      << ") has not the good number of component.");

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element = filter_elements->getSize();
    filter_elem_val = filter_elements->values;
  }

  out_uq.resize(nb_element);

  UInt * conn_val = mesh->getConnectivity(type).values;

  Real * shape_val = shapes[type]->values;
  Real * u_val     = in_u.values;
  Real * uq_val    = out_uq.values;

  UInt offset_uq   = out_uq.getNbComponent();

  Real * shape = shape_val;
  Real * u = static_cast<Real *>(calloc(nb_nodes_per_element * nb_degre_of_freedom,
					sizeof(Real)));

  for (UInt el = 0; el < nb_element; ++el) {
    UInt el_offset = el * nb_nodes_per_element;
    if(filter_elements != NULL) {
      shape     = shape_val + filter_elem_val[el] * size_of_shapes;
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
      shape += size_of_shapes;
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
				     const Vector<UInt> * filter_elements) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(shapes[type] != NULL,
		      "No shapes for the type " << type);

  UInt nb_nodes_per_element       = mesh->getConnectivity(type).getNbComponent();
  UInt size_of_shapes_derivatives = shapes_derivatives[type]->getNbComponent();
  UInt nb_quadrature_points       = getNbQuadraturePoints(type);
  UInt nb_element = mesh->getConnectivity(type).getSize();

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element = filter_elements->getSize();
    filter_elem_val = filter_elements->values;
  }

  AKANTU_DEBUG_ASSERT(in_u.getSize() == mesh->getNbNodes(),
		      "The vector in_u(" << in_u.getID()
		      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(in_u.getNbComponent() == nb_degre_of_freedom,
		      "The vector in_u(" << in_u.getID()
		      << ") has not the good number of component.");

  AKANTU_DEBUG_ASSERT(out_nablauq.getNbComponent() == nb_degre_of_freedom * nb_quadrature_points * spatial_dimension,
		      "The vector out_nablauq(" << out_nablauq.getID()
		      << ") has not the good number of component.");

  out_nablauq.resize(nb_element);

  UInt * conn_val = mesh->getConnectivity(type).values;

  Real * shaped_val  = shapes_derivatives[type]->values;
  Real * u_val       = in_u.values;
  Real * nablauq_val = out_nablauq.values;

  UInt offset_nablauq = nb_degre_of_freedom * spatial_dimension;
  UInt offset_shaped  = nb_nodes_per_element * spatial_dimension;

  Real * shaped  = shaped_val;
  Real * u       = static_cast<Real *>(calloc(nb_nodes_per_element * nb_degre_of_freedom,
					      sizeof(Real)));

  for (UInt el = 0; el < nb_element; ++el) {
    UInt el_offset = el * nb_nodes_per_element;
    if(filter_elements != NULL) {
      shaped    = shaped_val  + filter_elem_val[el] * size_of_shapes_derivatives;
      el_offset = filter_elem_val[el] * nb_nodes_per_element;
    }

    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      memcpy(u + n * nb_degre_of_freedom,
	     u_val + conn_val[el_offset + n] * nb_degre_of_freedom,
	     nb_degre_of_freedom * sizeof(Real));
    }

    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      /// \nabla(U) = U^t * dphi/dx
      Math::matrixt_matrix(nb_degre_of_freedom, spatial_dimension, nb_nodes_per_element,
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
		    const Vector<UInt> * filter_elements) {
  AKANTU_DEBUG_IN();


  UInt nb_element = filter_elements == NULL ? mesh->getNbElement(type) : filter_elements->getSize();
  //  UInt nb_nodes_per_element = mesh->getConnectivity(type).getNbComponent();
  UInt nb_quadrature_points = getNbQuadraturePoints(type);
  UInt size_of_jacobians    = jacobians[type]->getNbComponent();

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    filter_elem_val = filter_elements->values;
  }

  AKANTU_DEBUG_ASSERT(in_f.getSize() == mesh->getNbElement(type),
		      "The vector in_f(" << in_f.getID()
		      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(in_f.getNbComponent() == nb_degre_of_freedom * nb_quadrature_points,
		      "The vector in_f(" << in_f.getID()
		      << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(intf.getNbComponent() == nb_degre_of_freedom,
		      "The vector intf(" << intf.getID()
		      << ") has not the good number of component.");

  intf.resize(nb_element);

  Real * in_f_val = in_f.values;
  Real * intf_val = intf.values;
  Real * jac_val  = jacobians[type]->values;

  UInt offset_in_f = in_f.getNbComponent();
  UInt offset_intf = intf.getNbComponent();

  Real * jac      = jac_val;

  for (UInt el = 0; el < nb_element; ++el) {
    if(filter_elements != NULL) {
      jac      = jac_val  + filter_elem_val[el] * size_of_jacobians;
    }

    integrate(in_f_val, jac, intf_val, nb_degre_of_freedom, nb_quadrature_points);

    in_f_val += offset_in_f;
    intf_val += offset_intf;
    if(filter_elements == NULL) {
      jac      += size_of_jacobians;
    }
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
Real FEM::integrate(const Vector<Real> & in_f,
		    const ElementType & type,
		    const Vector<UInt> * filter_elements) {
  AKANTU_DEBUG_IN();

  UInt nb_element = filter_elements == NULL ? mesh->getNbElement(type) : filter_elements->getSize();
  //  UInt nb_nodes_per_element = mesh->getConnectivity(type).getNbComponent();
  UInt nb_quadrature_points = getNbQuadraturePoints(type);
  UInt size_of_jacobians    = jacobians[type]->getNbComponent();

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    filter_elem_val = filter_elements->values;
  }

  AKANTU_DEBUG_ASSERT(in_f.getSize() == mesh->getNbElement(type),
		      "The vector in_f(" << in_f.getID()
		      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(in_f.getNbComponent() == nb_quadrature_points,
		      "The vector in_f(" << in_f.getID()
		      << ") has not the good number of component.");

  Real intf = 0.;
  Real * in_f_val = in_f.values;
  Real * jac_val  = jacobians[type]->values;

  UInt offset_in_f = in_f.getNbComponent();

  Real * jac      = jac_val;

  for (UInt el = 0; el < nb_element; ++el) {
    if(filter_elements != NULL) {
      jac = jac_val  + filter_elem_val[el] * size_of_jacobians;
    }
    Real el_intf;
    integrate(in_f_val, jac, &el_intf, 1, nb_quadrature_points);
    intf += el_intf;

    in_f_val += offset_in_f;
    if(filter_elements == NULL) {
      jac += size_of_jacobians;
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
			 const Vector<UInt> * filter_elements,
			 Real scale_factor,
			 bool is_init_to_zero) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element = mesh->getNbNodesPerElement(type);
  UInt nb_element           = filter_elements == NULL ? mesh->getNbElement(type) : filter_elements->getSize();
  UInt nb_nodes = mesh->getNbNodes();

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    filter_elem_val = filter_elements->values;
  }

  AKANTU_DEBUG_ASSERT(elementary_vect.getSize() == nb_element,
		      "The vector elementary_vect(" << elementary_vect.getID()
		      << ") has not the good size.");

  AKANTU_DEBUG_ASSERT(elementary_vect.getNbComponent() == nb_degre_of_freedom * nb_nodes_per_element,
		      "The vector elementary_vect(" << elementary_vect.getID()
		      << ") has not the good number of component.");

  AKANTU_DEBUG_ASSERT(nodal_values.getNbComponent() == nb_degre_of_freedom,
		      "The vector nodal_values(" << nodal_values.getID()
		      << ") has not the good number of component.");

  nodal_values.resize(nb_nodes);

  UInt * conn_val = mesh->getConnectivity(type).values;
  Real * elementary_vect_val = elementary_vect.values;
  Real * nodal_values_val = nodal_values.values;

  if(!is_init_to_zero)
    memset(nodal_values_val, 0, nb_nodes*nb_degre_of_freedom*sizeof(Real));

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
  stream << space << " + spatial dimension : " << spatial_dimension << std::endl;

  stream << space << " + mesh [" << std::endl;
  mesh->printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << " + connectivity type information [" << std::endl;
  const Mesh::ConnectivityTypeList & type_list = mesh->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if (mesh->getSpatialDimension(*it) != spatial_dimension) continue;
    stream << space << AKANTU_INDENT << AKANTU_INDENT << " + " << *it <<" [" << std::endl;
    shapes            [*it]->printself(stream, indent + 3);
    shapes_derivatives[*it]->printself(stream, indent + 3);
    jacobians         [*it]->printself(stream, indent + 3);
    stream << space << AKANTU_INDENT << AKANTU_INDENT << "]" << std::endl;
  }
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << "]" << std::endl;
}


__END_AKANTU__
