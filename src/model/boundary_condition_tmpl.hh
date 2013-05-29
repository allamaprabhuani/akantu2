/**
 * @file   boundary_condition_inline_impl.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 *
 * @date   Wed Mar 18 11:30:00 2013
 *
 * @brief  XXX
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
template<typename ModelType>
void BoundaryCondition<ModelType>::initBC(ModelType & ref, Array<Real> & primal_, Array<Real> & dual_)
{
  model = &ref;
  primal = &primal_;
  dual = &dual_;
  model->initFEMBoundary();
}


/* -------------------------------------------------------------------------- */
/* Partial specialization for DIRICHLET functors */
template<typename ModelType>
template<typename FunctorType>
struct BoundaryCondition<ModelType>::TemplateFunctionWrapper<FunctorType, BC::Functor::_dirichlet> {
  static inline void applyBC(const FunctorType & func, const SubBoundary & boundary_ref, BoundaryCondition<ModelType> & bc_instance) {
    ModelType &         model          = *bc_instance.model;
    const Array<Real> & coords         = model.mesh.getNodes();
    UInt                dim            = model.mesh.getSpatialDimension();
    Array<Real>       & primal         = *bc_instance.primal;
    Array<bool>       & boundary_flags = model.getBoundary();

    Array<Real>::iterator<Vector<Real> > primal_iter = primal.begin(primal.getNbComponent());
    Array<Real>::const_iterator<Vector<Real> > coords_iter = coords.begin(dim);
    Array<bool>::iterator<Vector<bool> > flags_iter = boundary_flags.begin(boundary_flags.getNbComponent());

    for(SubBoundary::nodes_const_iterator nodes_it(boundary_ref.nodes_begin()); nodes_it!= boundary_ref.nodes_end(); ++nodes_it) {
      UInt n = *nodes_it;
      func(n, flags_iter[n], primal_iter[n], coords_iter[n]);
    }
  }
};


/* -------------------------------------------------------------------------- */
/* Partial specialization for NEUMANN functors */
template<typename ModelType>
template<typename FunctorType>
struct BoundaryCondition<ModelType>::TemplateFunctionWrapper<FunctorType, BC::Functor::_neumann> {
  static inline void applyBC(const FunctorType & func, const SubBoundary & boundary_ref, BoundaryCondition<ModelType> & bc_instance) {
    ModelType &                model        = *bc_instance.model;
    const Array<Real>        & nodes_coords = model.getMesh().getNodes();
    UInt                       dim          = model.getSpatialDimension();
    Array<Real>              & dual         = *bc_instance.dual;
    const FEM                & fem_boundary = model.getFEMBoundary();

    UInt nb_degree_of_freedom = dual.getNbComponent();

    Mesh::type_iterator type_it = model.mesh.firstType(dim-1);
    Mesh::type_iterator type_end = model.mesh.lastType(dim-1);

    // Loop over the boundary element types
    for(; type_it != type_end; ++type_it) {
      const Array<UInt> & element_ids = boundary_ref.getElements(*type_it);

      Array<UInt>::const_iterator<UInt> elem_iter = element_ids.begin();
      Array<UInt>::const_iterator<UInt> elem_iter_end = element_ids.end();

      UInt nb_quad_points = model.getFEMBoundary().getNbQuadraturePoints(*type_it);
      UInt nb_elements = element_ids.getSize();
      UInt nb_nodes_per_element = model.getMesh().getNbNodesPerElement(*type_it);

      Array<Real> * dual_before_integ = new Array<Real>(nb_elements * nb_quad_points, nb_degree_of_freedom);
      Array<Real> * quad_coords = new Array<Real>(nb_elements * nb_quad_points, dim);

      Array<Real>::iterator<Vector<Real> > dual_iter = dual_before_integ->begin(nb_degree_of_freedom);

      const Array<Real> & normals_on_quad = fem_boundary.getNormalsOnQuadPoints(*type_it);

      fem_boundary.interpolateOnQuadraturePoints(nodes_coords,
         *quad_coords, dim, *type_it, _not_ghost, element_ids);

      Array<Real>::const_iterator< Vector<Real> > normals_iter = normals_on_quad.begin(dim);
      Array<Real>::const_iterator< Vector<Real> > quad_coords_iter  = quad_coords->begin(dim);

      QuadraturePoint quad_point;
      quad_point.type = *type_it;

      // Loop over the elements in the SubBoundary
      for(; elem_iter != elem_iter_end; ++elem_iter) {
        UInt n = *elem_iter;
        quad_point.element = n;
        UInt offset = (n*nb_quad_points);
        for(UInt j(0); j<nb_quad_points; ++j) {

          quad_point.num_point = j;

          func(quad_point, *dual_iter, *quad_coords_iter, normals_iter[offset]);

          ++dual_iter;
          ++quad_coords_iter;
          ++offset;
        }
      }
      delete quad_coords;

      Array<Real>::iterator<Matrix<Real> > dual_iter_mat = dual_before_integ->begin(nb_degree_of_freedom,1);
      elem_iter = element_ids.begin();
      Array<Real>::const_iterator<Matrix<Real> > shapes_iter_begin = fem_boundary.getShapes(*type_it).begin(1, nb_nodes_per_element);

      Array<Real> * dual_by_shapes = new Array<Real>(nb_elements*nb_quad_points, nb_degree_of_freedom*nb_nodes_per_element);
      Array<Real>::iterator<Matrix<Real> > dual_by_shapes_iter = dual_by_shapes->begin(nb_degree_of_freedom, nb_nodes_per_element);

      for(; elem_iter != elem_iter_end; ++elem_iter) {
        Array<Real>::const_iterator<Matrix<Real> > shapes_iter = shapes_iter_begin + *elem_iter*nb_quad_points;

        for(UInt j(0); j<nb_quad_points; ++j, ++dual_iter_mat, ++dual_by_shapes_iter, ++shapes_iter) {
          dual_by_shapes_iter->mul<false, false>(*dual_iter_mat, *shapes_iter);
        }
      }
      delete dual_before_integ;
      Array<Real> * dual_by_shapes_integ = new Array<Real>(nb_elements, nb_degree_of_freedom*nb_nodes_per_element);

      fem_boundary.integrate(*dual_by_shapes, *dual_by_shapes_integ, nb_degree_of_freedom*nb_nodes_per_element, *type_it, _not_ghost, element_ids);
      delete dual_by_shapes;

      // assemble the result into force vector
      fem_boundary.assembleArray(*dual_by_shapes_integ, dual, model.getDOFSynchronizer().getLocalDOFEquationNumbers(), nb_degree_of_freedom, *type_it, _not_ghost, element_ids);
      delete dual_by_shapes_integ;
    }
  }
};

/* -------------------------------------------------------------------------- */
template<typename ModelType>
template<typename FunctorType>
inline void BoundaryCondition<ModelType>::applyBC(const FunctorType & func) {
  Boundary::const_iterator bit = model->getMesh().getBoundary().begin();
  Boundary::const_iterator bend = model->getMesh().getBoundary().end();
  for(; bit != bend; ++bit) {
    TemplateFunctionWrapper<FunctorType>::applyBC(func, *bit, *this);
  }
}

/* -------------------------------------------------------------------------- */
template<typename ModelType>
template<typename FunctorType>
inline void BoundaryCondition<ModelType>::applyBC(const FunctorType & func, const std::string & boundary_name) {
  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  const SubBoundary * boundary_ptr = NULL;
  try {
    boundary_ptr = &(model->getMesh().getSubBoundary(boundary_name));
    TemplateFunctionWrapper<FunctorType>::applyBC(func, *boundary_ptr, *this);
  } catch(akantu::debug::Exception e) {
    if(psize == 1) {
      AKANTU_DEBUG_ERROR("Error applying a boundary condition onto \"" << boundary_name << "\" in  a serial program. This should not occur!");
    }
  }
}

/* -------------------------------------------------------------------------- */
template<typename ModelType>
template<typename FunctorType>
inline void BoundaryCondition<ModelType>::applyBC(const FunctorType & func, const SubBoundary & boundary_ref) {
  TemplateFunctionWrapper<FunctorType>::applyBC(func, boundary_ref, *this);
}
