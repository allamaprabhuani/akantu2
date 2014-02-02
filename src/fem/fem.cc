/**
 * @file   fem.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Jul 20 23:40:43 2010
 *
 * @brief  Implementation of the FEM class
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
#include "fem.hh"
#include "mesh.hh"
#include "element_class.hh"
#include "static_communicator.hh"
#include "aka_math.hh"
#include "dof_synchronizer.hh"


/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
FEM::FEM(Mesh & mesh, UInt element_dimension, ID id, MemoryID memory_id) :
  Memory(id, memory_id), mesh(mesh), normals_on_quad_points("normals_on_quad_points", id) {
  AKANTU_DEBUG_IN();
  this->element_dimension = (element_dimension != _all_dimensions) ?
    element_dimension : mesh.getSpatialDimension();

  init();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FEM::init() {

}

/* -------------------------------------------------------------------------- */
FEM::~FEM() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void FEM::assembleArray(const Array<Real> & elementary_vect,
			 Array<Real> & nodal_values,
			 const Array<Int> & equation_number,
			 UInt nb_degree_of_freedom,
			 const ElementType & type,
			 const GhostType & ghost_type,
			 const Array<UInt> & filter_elements,
			 Real scale_factor) const {
  AKANTU_DEBUG_IN();

  UInt nb_element;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  Array<UInt>::const_iterator< Vector<UInt> > conn_it;

  Array<UInt> * filtered_connectivity = NULL;
  if(filter_elements != empty_filter) {
    nb_element = filter_elements.getSize();
    filtered_connectivity = new Array<UInt>(0, nb_nodes_per_element);
    FEM::filterElementalData(mesh,
                             mesh.getConnectivity(type, ghost_type),
                             *filtered_connectivity,
                             type, ghost_type,
                             filter_elements);
    const Array<UInt> & cfiltered = *filtered_connectivity; // \todo temporary patch
    conn_it = cfiltered.begin(nb_nodes_per_element);
  } else {
    nb_element = mesh.getNbElement(type, ghost_type);
    conn_it = mesh.getConnectivity(type, ghost_type).begin(nb_nodes_per_element);
  }

  AKANTU_DEBUG_ASSERT(elementary_vect.getSize() == nb_element,
		      "The vector elementary_vect(" << elementary_vect.getID()
		      << ") has not the good size.");

  AKANTU_DEBUG_ASSERT(elementary_vect.getNbComponent()
		      == nb_degree_of_freedom*nb_nodes_per_element,
		      "The vector elementary_vect(" << elementary_vect.getID()
		      << ") has not the good number of component."
		      << "(" << elementary_vect.getNbComponent()
		      << " != " << nb_degree_of_freedom*nb_nodes_per_element << ")");

  AKANTU_DEBUG_ASSERT(nodal_values.getNbComponent() == nb_degree_of_freedom,
		      "The vector nodal_values(" << nodal_values.getID()
		      << ") has not the good number of component."
		      << "(" << nodal_values.getNbComponent()
		      << " != " << nb_degree_of_freedom << ")");


  nodal_values.resize(mesh.getNbNodes());
  Real * nodal_it  = nodal_values.storage();
  Array<Real>::const_iterator< Matrix<Real> > elem_it  = elementary_vect.begin(nb_degree_of_freedom,
                                                                               nb_nodes_per_element);

  for (UInt el = 0; el < nb_element; ++el, ++elem_it, ++conn_it) {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt node = (*conn_it)(n);
      UInt offset_node = node * nb_degree_of_freedom;

      const Vector<Real> & elem_data = (*elem_it)(n);
      for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
        nodal_it[equation_number(offset_node + d)]
          += scale_factor * elem_data(d);
      }
    }
  }

  delete filtered_connectivity;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FEM::assembleMatrix(const Array<Real> & elementary_mat,
			 SparseMatrix & matrix,
			 UInt nb_degree_of_freedom,
			 const ElementType & type,
			 const GhostType & ghost_type,
			 const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();


  UInt nb_element;
  if(ghost_type == _not_ghost) {
    nb_element  = mesh.getNbElement(type);
  } else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  if(filter_elements != empty_filter) {
    nb_element      = filter_elements.getSize();
  }

  AKANTU_DEBUG_ASSERT(elementary_mat.getSize() == nb_element,
		      "The vector elementary_mat(" << elementary_mat.getID()
		      << ") has not the good size.");

  AKANTU_DEBUG_ASSERT(elementary_mat.getNbComponent()
		      == nb_degree_of_freedom * nb_nodes_per_element * nb_degree_of_freedom * nb_nodes_per_element,
		      "The vector elementary_mat(" << elementary_mat.getID()
		      << ") has not the good number of component.");

  Real * elementary_mat_val = elementary_mat.values;
  UInt offset_elementary_mat = elementary_mat.getNbComponent();
  UInt * connectivity_val = mesh.getConnectivity(type, ghost_type).values;

  UInt size_mat = nb_nodes_per_element * nb_degree_of_freedom;
  UInt size = mesh.getNbGlobalNodes() * nb_degree_of_freedom;

  Int * eq_nb_val = matrix.getDOFSynchronizer().getGlobalDOFEquationNumbers().values;
  Int * local_eq_nb_val = new Int[size_mat];

  for (UInt e = 0; e < nb_element; ++e) {
    UInt el = e;
    if(filter_elements != empty_filter) el = filter_elements(e);

    Int * tmp_local_eq_nb_val = local_eq_nb_val;
    UInt * conn_val = connectivity_val + el * nb_nodes_per_element;
    for (UInt i = 0; i < nb_nodes_per_element; ++i) {
      UInt n = conn_val[i];
      for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
        *tmp_local_eq_nb_val++ = eq_nb_val[n * nb_degree_of_freedom + d];
      }
      // memcpy(tmp_local_eq_nb_val, eq_nb_val + n * nb_degree_of_freedom, nb_degree_of_freedom * sizeof(Int));
      // tmp_local_eq_nb_val += nb_degree_of_freedom;
    }

    for (UInt i = 0; i < size_mat; ++i) {
      UInt c_irn = local_eq_nb_val[i];
      if(c_irn < size) {
        UInt j_start = (matrix.getSparseMatrixType() == _symmetric) ? i : 0;
        for (UInt j = j_start; j < size_mat; ++j) {
          UInt c_jcn = local_eq_nb_val[j];
          if(c_jcn < size) {
            matrix(c_irn, c_jcn) += elementary_mat_val[j * size_mat + i];
          }
        }
      }
    }
    elementary_mat_val += offset_elementary_mat;
  }

  delete [] local_eq_nb_val;

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
  mesh.printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;


  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */



__END_AKANTU__
