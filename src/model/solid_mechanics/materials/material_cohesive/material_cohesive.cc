/**
 * @file   material_cohesive.cc
 *
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Wed Feb 22 16:31:20 2012
 *
 * @brief  Specialization of the material class for cohesive elements
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
#include "material_cohesive.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "sparse_matrix.hh"
#include "dof_synchronizer.hh"
#include "aka_random_generator.hh"
#include "shape_cohesive.hh"


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialCohesive::MaterialCohesive(SolidMechanicsModel & model, const ID & id) :
  Material(model,id),
  reversible_energy("reversible_energy", id),
  total_energy("total_energy", id),
  tractions_old("tractions (old)",id),
  opening_old("opening (old)",id),
  tractions("tractions",id),
  opening("opening",id),
  contact_tractions("contact_tractions",id),
  contact_opening("contact_opening",id),
  delta_max("delta max",id),
  damage("damage", id),
  facet_filter("facet_filter", id),
  sigma_limit("sigma_limit", id) {

  AKANTU_DEBUG_IN();

  this->model = dynamic_cast<SolidMechanicsModelCohesive*>(&model);

  fem_cohesive = &(model.getFEMClass<MyFEMCohesiveType>("CohesiveFEM"));

  this->registerParam("sigma_c", sigma_c, 0.,
		      ParamAccessType(_pat_parsable | _pat_readable),
		      "Critical stress");

  this->registerParam("sigma_c_generator", random_generator,
		      _pat_parsable, "Random generator for sigma_c");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MaterialCohesive::~MaterialCohesive() {
  AKANTU_DEBUG_IN();

  delete random_generator;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::initMaterial() {
  AKANTU_DEBUG_IN();

  Material::initMaterial();

  initInternalArray(reversible_energy,                 1, false, _ek_cohesive);
  initInternalArray(     total_energy,                 1, false, _ek_cohesive);
  initInternalArray(    tractions_old, spatial_dimension, false, _ek_cohesive);
  initInternalArray(        tractions, spatial_dimension, false, _ek_cohesive);
  initInternalArray(      opening_old, spatial_dimension, false, _ek_cohesive);
  initInternalArray(contact_tractions, spatial_dimension, false, _ek_cohesive);
  initInternalArray(  contact_opening, spatial_dimension, false, _ek_cohesive);
  initInternalArray(          opening, spatial_dimension, false, _ek_cohesive);
  initInternalArray(        delta_max,                 1, false, _ek_cohesive);
  initInternalArray(           damage,                 1, false, _ek_cohesive);
  initInternalArray(   element_filter,                 1, false, _ek_cohesive);

  AKANTU_DEBUG_OUT();
}

void MaterialCohesive::initInsertionArrays(const Mesh & mesh_facets) {
  AKANTU_DEBUG_IN();

  mesh_facets.initByElementTypeArray(facet_filter, 1, spatial_dimension - 1);
  mesh_facets.initByElementTypeArray( sigma_limit, 1, spatial_dimension - 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::resizeCohesiveArrays() {
  AKANTU_DEBUG_IN();

  resizeInternalArray(reversible_energy, _ek_cohesive);
  resizeInternalArray(total_energy     , _ek_cohesive);
  resizeInternalArray(tractions_old    , _ek_cohesive);
  resizeInternalArray(tractions        , _ek_cohesive);
  resizeInternalArray(opening_old      , _ek_cohesive);
  resizeInternalArray(opening          , _ek_cohesive);
  resizeInternalArray(contact_tractions, _ek_cohesive);
  resizeInternalArray(contact_opening  , _ek_cohesive);
  resizeInternalArray(delta_max        , _ek_cohesive);
  resizeInternalArray(damage           , _ek_cohesive);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void MaterialCohesive::generateRandomDistribution(const Mesh & mesh_facets) {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it   = mesh_facets.firstType(spatial_dimension - 1);
  Mesh::type_iterator last = mesh_facets.lastType(spatial_dimension - 1);

  for (; it != last; ++it) {
    ElementType type_facet = *it;
    Array<UInt> & f_filter = facet_filter(type_facet);
    UInt nb_facet = f_filter.getSize();

    if (nb_facet > 0) {
      Array<Real> & sigma_lim = sigma_limit(type_facet);
      sigma_lim.resize(nb_facet);

      if (random_generator)
	random_generator->generate(sigma_lim);
      else
      	sigma_lim.set(sigma_c);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::assembleResidual(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model->getSpatialDimension();

  Array<Real> & residual = const_cast<Array<Real> &>(model->getResidual());

  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension,
					  ghost_type, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension,
						ghost_type, _ek_cohesive);
  for(; it != last_type; ++it) {
    const Array<Real> & shapes = fem_cohesive->getShapes(*it, ghost_type);

    Array<UInt> & elem_filter = element_filter(*it, ghost_type);
    Array<Real> & traction = tractions(*it, ghost_type);

    UInt size_of_shapes       = shapes.getNbComponent();
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
    UInt nb_quadrature_points = fem_cohesive->getNbQuadraturePoints(*it, ghost_type);

    UInt nb_element = elem_filter.getSize();
    if (nb_element == 0) continue;

    /// compute @f$t_i N_a@f$

    Real * shapes_val       = shapes.storage();
    UInt * elem_filter_val  = elem_filter.storage();

    Array<Real> * shapes_filtered = new Array<Real>(nb_element*nb_quadrature_points,
						    size_of_shapes,
						    "filtered shapes");

    Real * shapes_filtered_val = shapes_filtered->values;

    for (UInt el = 0; el < nb_element; ++el) {
      shapes_val = shapes.storage() + elem_filter_val[el] *
	size_of_shapes * nb_quadrature_points;
      memcpy(shapes_filtered_val, shapes_val,
	     size_of_shapes * nb_quadrature_points * sizeof(Real));
      shapes_filtered_val += size_of_shapes * nb_quadrature_points;
    }

    shapes_filtered_val = shapes_filtered->values;

    // multiply traction by shapes

    Array<Real> * traction_cpy = new Array<Real>(traction);

    /// add contact tractions
    Array<Real> & contact_traction = contact_tractions(*it, ghost_type);
    (*traction_cpy) += contact_traction;

    traction_cpy->extendComponentsInterlaced(size_of_shapes, spatial_dimension);

    Real * traction_cpy_val = traction_cpy->storage();

    for (UInt el = 0; el < nb_element; ++el) {
      for (UInt q = 0; q < nb_quadrature_points; ++q) {
	for (UInt n = 0; n < size_of_shapes; ++n,++shapes_filtered_val) {
	  for (UInt i = 0; i < spatial_dimension; ++i) {
	    *traction_cpy_val++ *= *shapes_filtered_val;
	  }
	}
      }
    }

    delete shapes_filtered;

    /**
     * compute @f$\int t \cdot N\, dS@f$ by  @f$ \sum_q \mathbf{N}^t
     * \mathbf{t}_q \overline w_q J_q@f$
     */
    Array<Real> * int_t_N = new Array<Real>(nb_element,
					    spatial_dimension*size_of_shapes,
					    "int_t_N");

    fem_cohesive->integrate(*traction_cpy, *int_t_N,
			    spatial_dimension*size_of_shapes,
			    *it, ghost_type,
			    elem_filter);

    delete traction_cpy;

    int_t_N->extendComponentsInterlaced(2, int_t_N->getNbComponent() );

    Real * int_t_N_val = int_t_N->storage();
    for (UInt el = 0; el < nb_element; ++el) {
      for (UInt n = 0; n < size_of_shapes*spatial_dimension; ++n)
	int_t_N_val[n] *= -1.;
      int_t_N_val += nb_nodes_per_element*spatial_dimension;
    }

    /// assemble
    model->getFEMBoundary().assembleArray(*int_t_N, residual,
					  model->getDOFSynchronizer().getLocalDOFEquationNumbers(),
					  residual.getNbComponent(),
					  *it, ghost_type, elem_filter, 1);

    delete int_t_N;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::assembleStiffnessMatrix(GhostType ghost_type) {

  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model->getSpatialDimension();

  SparseMatrix & K = const_cast<SparseMatrix &>(model->getStiffnessMatrix());

  Mesh & mesh = fem_cohesive->getMesh();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension,
					  ghost_type, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension,
						ghost_type, _ek_cohesive);

  for(; it != last_type; ++it) {
    UInt nb_quadrature_points = fem_cohesive->getNbQuadraturePoints(*it, ghost_type);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);

    const Array<Real> & shapes = fem_cohesive->getShapes(*it, ghost_type);
    Array<UInt> & elem_filter = element_filter(*it, ghost_type);

    UInt nb_element = elem_filter.getSize();
    UInt size_of_shapes       = shapes.getNbComponent();

    Array<Real> * shapes_filtered =
      new Array<Real>(nb_element*nb_quadrature_points,
		      size_of_shapes, "filtered shapes");

    Real * shapes_val       = shapes.storage();
    Real * shapes_filtered_val = shapes_filtered->values;
    UInt * elem_filter_val  = elem_filter.storage();

    for (UInt el = 0; el < nb_element; ++el) {
      shapes_val = shapes.storage() + elem_filter_val[el] *
    	size_of_shapes * nb_quadrature_points;
      memcpy(shapes_filtered_val, shapes_val,
    	     size_of_shapes * nb_quadrature_points * sizeof(Real));
      shapes_filtered_val += size_of_shapes * nb_quadrature_points;
    }

    /**
     * compute A matrix @f$ \mathbf{A} = \left[\begin{array}{c c c c c c c c c c c c}
     * 1 & 0 & 0 & 0 & 0 & 0 & -1 &  0 &  0 &  0 &  0 &  0 \\
     * 0 & 1 & 0 & 0 & 0 & 0 &  0 & -1 &  0 &  0 &  0 &  0 \\
     * 0 & 0 & 1 & 0 & 0 & 0 &  0 &  0 & -1 &  0 &  0 &  0 \\
     * 0 & 0 & 0 & 1 & 0 & 0 &  0 &  0 &  0 & -1 &  0 &  0 \\
     * 0 & 0 & 0 & 0 & 1 & 0 &  0 &  0 &  0 &  0 & -1 &  0 \\
     * 0 & 0 & 0 & 0 & 0 & 1 &  0 &  0 &  0 &  0 &  0 & -1
     * \end{array} \right]@f$
     **/

    // UInt size_of_A =  spatial_dimension*size_of_shapes*spatial_dimension*nb_nodes_per_element;
    // Real * A = new Real[size_of_A];
    // memset(A, 0, size_of_A*sizeof(Real));
    Matrix<Real> A(spatial_dimension*size_of_shapes,
		   spatial_dimension*nb_nodes_per_element);

    for ( UInt i = 0; i < spatial_dimension*size_of_shapes; ++i) {
      A(i, i) = 1;
      A(i, i + spatial_dimension*size_of_shapes) = -1;
    }

    /// compute traction
    computeTraction(ghost_type);

    /// get the tangent matrix @f$\frac{\partial{(t/\delta)}}{\partial{\delta}} @f$
    Array<Real> * tangent_stiffness_matrix =
      new Array<Real>(nb_element * nb_quadrature_points,
		      spatial_dimension * spatial_dimension,
		      "tangent_stiffness_matrix");

    //    Array<Real> * normal = new Array<Real>(nb_element * nb_quadrature_points, spatial_dimension, "normal");
    Array<Real> normal(nb_quadrature_points, spatial_dimension, "normal");


    computeNormal(model->getCurrentPosition(), normal, *it, ghost_type);

    tangent_stiffness_matrix->clear();

    computeTangentTraction(*it, *tangent_stiffness_matrix, normal, ghost_type);

    // delete normal;

    UInt size_at_nt_d_n_a = spatial_dimension*nb_nodes_per_element*spatial_dimension*nb_nodes_per_element;
    Array<Real> * at_nt_d_n_a = new Array<Real> (nb_element*nb_quadrature_points,
						   size_at_nt_d_n_a,
						   "A^t*N^t*D*N*A");

    Array<Real>::iterator<Vector<Real> > shapes_filt_it =
      shapes_filtered->begin(size_of_shapes);

    Array<Real>::iterator< Matrix<Real> > D_it =
      tangent_stiffness_matrix->begin(spatial_dimension, spatial_dimension);

    Array<Real>::iterator< Matrix<Real> > At_Nt_D_N_A_it =
      at_nt_d_n_a->begin(spatial_dimension * nb_nodes_per_element,
			 spatial_dimension * nb_nodes_per_element);

    Array<Real>::iterator< Matrix<Real> > At_Nt_D_N_A_end =
      at_nt_d_n_a->end(spatial_dimension * nb_nodes_per_element,
		       spatial_dimension * nb_nodes_per_element);

    Matrix<Real> N    (spatial_dimension, spatial_dimension * size_of_shapes);
    Matrix<Real> N_A  (spatial_dimension, spatial_dimension * nb_nodes_per_element);
    Matrix<Real> D_N_A(spatial_dimension, spatial_dimension * nb_nodes_per_element);

    for(; At_Nt_D_N_A_it != At_Nt_D_N_A_end; ++At_Nt_D_N_A_it, ++D_it, ++shapes_filt_it) {
      N.clear();
      /**
       * store  the   shapes  in  voigt   notations  matrix  @f$\mathbf{N}  =
       * \begin{array}{cccccc} N_0(\xi) & 0 & N_1(\xi)  &0 & N_2(\xi) & 0 \\
       * 0 & * N_0(\xi)& 0 &N_1(\xi)& 0 & N_2(\xi) \end{array} @f$
       **/
      for (UInt i = 0; i < spatial_dimension ; ++i)
	for (UInt n = 0; n < size_of_shapes; ++n)
	  N(i, i + spatial_dimension * n) = (*shapes_filt_it)(n);

      /**
       * compute stiffness matrix  @f$   \mathbf{K}    =    \delta    \mathbf{U}^T
       * \int_{\Gamma_c}    {\mathbf{P}^t \frac{\partial{\mathbf{t}}} {\partial{\delta}}
       * \mathbf{P} d\Gamma \Delta \mathbf{U}}  @f$
       **/
      N_A.mul<false, false>(N, A);
      D_N_A.mul<false, false>(*D_it, N_A);
      (*At_Nt_D_N_A_it).mul<true, false>(D_N_A, N_A);
    }

    delete tangent_stiffness_matrix;
    delete shapes_filtered;

    Array<Real> * K_e = new Array<Real>(nb_element, size_at_nt_d_n_a,
					  "K_e");

    fem_cohesive->integrate(*at_nt_d_n_a, *K_e,
			    size_at_nt_d_n_a,
			    *it, ghost_type,
			    elem_filter);

    delete at_nt_d_n_a;

    model->getFEM().assembleMatrix(*K_e, K, spatial_dimension,
				   *it, ghost_type, elem_filter);
    delete K_e;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- *
 * Compute traction from displacements
 *
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void MaterialCohesive::computeTraction(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model->getSpatialDimension();

  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension,
					  ghost_type, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension,
						ghost_type, _ek_cohesive);

  for(; it != last_type; ++it) {
    Array<UInt> & elem_filter = element_filter(*it, ghost_type);

    UInt nb_element = elem_filter.getSize();
    if (nb_element == 0) continue;
    UInt nb_quadrature_points =
      nb_element*fem_cohesive->getNbQuadraturePoints(*it, ghost_type);

    Array<Real> normal(nb_quadrature_points, spatial_dimension, "normal");

    /// compute normals @f$\mathbf{n}@f$
    computeNormal(model->getCurrentPosition(), normal, *it, ghost_type);

    /// compute openings @f$\mathbf{\delta}@f$
    computeOpening(model->getDisplacement(), opening(*it, ghost_type), *it, ghost_type);

    /// compute traction @f$\mathbf{t}@f$
    computeTraction(normal, *it, ghost_type);

  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void MaterialCohesive::computeNormal(const Array<Real> & position,
				     Array<Real> & normal,
				     ElementType type,
				     GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if (type == _cohesive_1d_2)
    fem_cohesive->computeNormalsOnControlPoints(position,
						normal,
						type, ghost_type);
  else {
#define COMPUTE_NORMAL(type)						\
    fem_cohesive->getShapeFunctions().					\
      computeNormalsOnControlPoints<type, CohesiveReduceFunctionMean>(position, \
								      normal, \
								      ghost_type, \
								      element_filter(type, ghost_type));

    AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_NORMAL);
#undef COMPUTE_NORMAL
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void MaterialCohesive::computeOpening(const Array<Real> & displacement,
				      Array<Real> & opening,
				      ElementType type,
				      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

#define COMPUTE_OPENING(type)						\
  fem_cohesive->getShapeFunctions().					\
    interpolateOnControlPoints<type, CohesiveReduceFunctionOpening>(displacement, \
								    opening, \
								    spatial_dimension, \
								    ghost_type, \
								    element_filter(type, ghost_type));

  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_OPENING);
#undef COMPUTE_OPENING

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::computeEnergies() {
  AKANTU_DEBUG_IN();

  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it =
    mesh.firstType(spatial_dimension, _not_ghost, _ek_cohesive);
  Mesh::type_iterator last_type =
    mesh.lastType(spatial_dimension, _not_ghost, _ek_cohesive);

  Real * memory_space = new Real[2*spatial_dimension];
  Vector<Real> b(memory_space, spatial_dimension);
  Vector<Real> h(memory_space + spatial_dimension, spatial_dimension);

  for(; it != last_type; ++it) {
    Array<Real>::iterator<Real> erev =
      reversible_energy(*it, _not_ghost).begin();
    Array<Real>::iterator<Real> etot =
      total_energy(*it, _not_ghost).begin();
    Array<Real>::iterator< Vector<Real> > traction_it =
      tractions(*it, _not_ghost).begin(spatial_dimension);
    Array<Real>::iterator< Vector<Real> > traction_old_it =
      tractions_old(*it, _not_ghost).begin(spatial_dimension);
    Array<Real>::iterator< Vector<Real> > opening_it =
      opening(*it, _not_ghost).begin(spatial_dimension);
    Array<Real>::iterator< Vector<Real> > opening_old_it =
      opening_old(*it, _not_ghost).begin(spatial_dimension);

    Array<Real>::iterator< Vector<Real> >traction_end =
      tractions(*it, _not_ghost).end(spatial_dimension);

    /// loop on each quadrature point
    for (; traction_it != traction_end;
	 ++traction_it, ++traction_old_it,
	   ++opening_it, ++opening_old_it,
	   ++erev, ++etot) {

      /// trapezoidal integration
      b  = *opening_it;
      b -= *opening_old_it;

      h  = *traction_old_it;
      h += *traction_it;

      *etot += .5 * b.dot(h);
      *erev  = .5 * traction_it->dot(*opening_it);
    }
  }

  delete [] memory_space;

  /// update old values
  it = mesh.firstType(spatial_dimension, _not_ghost, _ek_cohesive);
  GhostType ghost_type = _not_ghost;
  for(; it != last_type; ++it) {
    tractions_old(*it, ghost_type).copy(tractions(*it, ghost_type));
    opening_old(*it, ghost_type).copy(opening(*it, ghost_type));
  }


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getReversibleEnergy() {
  AKANTU_DEBUG_IN();
  Real erev = 0.;

  /// integrate reversible energy for each type of elements
  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension,
					  _not_ghost, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension,
						_not_ghost, _ek_cohesive);

  for(; it != last_type; ++it) {
    erev += fem_cohesive->integrate(reversible_energy(*it, _not_ghost), *it,
				    _not_ghost, element_filter(*it, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return erev;
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getDissipatedEnergy() {
  AKANTU_DEBUG_IN();
  Real edis = 0.;

  /// integrate dissipated energy for each type of elements
  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension,
					  _not_ghost, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension,
						_not_ghost, _ek_cohesive);

  for(; it != last_type; ++it) {
    Array<Real> dissipated_energy(total_energy(*it, _not_ghost));
    dissipated_energy -= reversible_energy(*it, _not_ghost);
    edis += fem_cohesive->integrate(dissipated_energy, *it,
				    _not_ghost, element_filter(*it, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return edis;
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getContactEnergy() {
  AKANTU_DEBUG_IN();
  Real econ = 0.;

  /// integrate contact energy for each type of elements
  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension,
					  _not_ghost, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension,
						_not_ghost, _ek_cohesive);

  for(; it != last_type; ++it) {
    Array<UInt> & el_filter = element_filter(*it, _not_ghost);
    UInt nb_element = el_filter.getSize();
    Array<Real> contact_energy(nb_element);

    Array<Real>::iterator< Vector<Real> > contact_traction_it =
      contact_tractions(*it, _not_ghost).begin(spatial_dimension);
    Array<Real>::iterator< Vector<Real> > contact_opening_it =
      contact_opening(*it, _not_ghost).begin(spatial_dimension);

    /// loop on each quadrature point
    for (UInt el = 0; el < nb_element;
	 ++contact_traction_it, ++contact_opening_it, ++el) {

      contact_energy(el) = .5 * contact_traction_it->dot(*contact_opening_it);
    }

    econ += fem_cohesive->integrate(contact_energy, *it,
				    _not_ghost, el_filter);
  }

  AKANTU_DEBUG_OUT();
  return econ;
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getEnergy(std::string type) {
  AKANTU_DEBUG_IN();

  if (type == "reversible") return getReversibleEnergy();
  else if (type == "dissipated") return getDissipatedEnergy();
  else if (type == "cohesive contact") return getContactEnergy();

  AKANTU_DEBUG_OUT();
  return 0.;
}


__END_AKANTU__
