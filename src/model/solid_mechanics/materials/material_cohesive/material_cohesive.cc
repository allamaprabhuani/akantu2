/**
 * @file   material_cohesive.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @date   Tue Feb  7 18:24:52 2012
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
  delta_max("delta max",id),
  damage("damage", id) {

  AKANTU_DEBUG_IN();

  this->model = dynamic_cast<SolidMechanicsModelCohesive*>(&model);

  initInternalVector(reversible_energy,                 1, false, _ek_cohesive);
  initInternalVector(     total_energy,                 1, false, _ek_cohesive);
  initInternalVector(    tractions_old, spatial_dimension, false, _ek_cohesive);
  initInternalVector(        tractions, spatial_dimension, false, _ek_cohesive);
  initInternalVector(      opening_old, spatial_dimension, false, _ek_cohesive);
  initInternalVector(          opening, spatial_dimension, false, _ek_cohesive);
  initInternalVector(        delta_max,                 1, false, _ek_cohesive);
  initInternalVector(           damage,                 1, false, _ek_cohesive);

  this->registerParam("sigma_c", sigma_c, 0. , _pat_parsable, "Critical stress");
  this->registerParam("rand_factor", rand, 0. , _pat_parsable, "Randomness factor");
  this->registerParam("distribution", distribution, std::string("uniform"),  _pat_parsable, "Distribution type");
  this->registerParam("lambda", lambda, 0. , _pat_parsable, "Weibull modulus");
  this->registerParam("m", m_scale, 1. , _pat_parsable, "Scale parameter");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MaterialCohesive::~MaterialCohesive() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::initMaterial() {
  AKANTU_DEBUG_IN();

  Material::initMaterial();

  fem_cohesive = &(model->getFEMClass<MyFEMCohesiveType>("CohesiveFEM"));

  Mesh & mesh = fem_cohesive->getMesh();

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    Mesh::type_iterator it = mesh.firstType(spatial_dimension, gt, _ek_cohesive);
    Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, gt, _ek_cohesive);

    for(; it != last_type; ++it) {
      ElementType type = *it;
      element_filter.alloc(0, 1, type);
    }
  }

  resizeCohesiveVectors();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::resizeCohesiveVectors() {
  resizeInternalVector(reversible_energy, _ek_cohesive);
  resizeInternalVector(total_energy, _ek_cohesive);
  resizeInternalVector(tractions_old, _ek_cohesive);
  resizeInternalVector(tractions, _ek_cohesive);
  resizeInternalVector(opening_old, _ek_cohesive);
  resizeInternalVector(opening, _ek_cohesive);
  resizeInternalVector(delta_max, _ek_cohesive);
  resizeInternalVector(damage, _ek_cohesive);
}

/* -------------------------------------------------------------------------- */

void MaterialCohesive::generateRandomDistribution(Vector<Real> & sigma_lim) {
  AKANTU_DEBUG_IN();

  std::srand(time(NULL));
  UInt nb_facet = sigma_lim.getSize();

  if (distribution == "uniform") {
    for (UInt i = 0; i < nb_facet; ++i)
      sigma_lim(i) = sigma_c * (1 + std::rand()/(Real)RAND_MAX * rand);
  }
  else if (distribution == "weibull") {
    Real exponent = 1./m_scale;
    for (UInt i = 0; i < nb_facet; ++i)
      sigma_lim(i) = sigma_c + lambda * std::pow(-1.* std::log(std::rand()/(Real)RAND_MAX), exponent);
  }
  else {
    AKANTU_DEBUG_ERROR("Unknown random distribution type for sigma_c");
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::checkInsertion(const Vector<Real> & facet_stress,
				      Vector<UInt> & facet_insertion) {
  AKANTU_DEBUG_IN();

  Vector<bool> & facets_check = model->getFacetsCheck();
  ElementType type_facet = model->getFacetType();

  UInt nb_quad_facet = model->getFEM("FacetsFEM").getNbQuadraturePoints(type_facet);
  UInt nb_facet = facets_check.getSize();

  Vector<Real> stress_check(nb_facet, nb_quad_facet);
  stress_check.clear();

  computeStressNorms(facet_stress, stress_check);

  bool * facet_check_it = facets_check.storage();
  Vector<Real>::iterator<types::RVector> stress_check_it =
    stress_check.begin(nb_quad_facet);

  Real * sigma_limit_it = model->getSigmaLimit().storage();

  for (UInt f = 0; f < nb_facet;
       ++f, ++facet_check_it, ++stress_check_it, ++sigma_limit_it) {
    if (*facet_check_it == true) {

      for (UInt q = 0; q < nb_quad_facet; ++q) {
  	if ((*stress_check_it)(q) > *sigma_limit_it) {
  	  facet_insertion.push_back(f);
  	  for (UInt qs = 0; qs < nb_quad_facet; ++qs)
  	    sigma_insertion.push_back((*stress_check_it)(qs));
  	  *facet_check_it = false;
  	  break;
  	}
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute  the  residual  by  assembling  @f$\int_{e}  t_e  N_e dS @f$
 *
 * @param[in] displacements nodes displacements
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void MaterialCohesive::updateResidual(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// compute traction
  computeTraction(ghost_type);

  /// update and assemble residual
  assembleResidual(ghost_type);

  /// compute energies
  computeEnergies();

  /// update old values
  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type, _ek_cohesive);
  for(; it != last_type; ++it) {
    tractions_old(*it, ghost_type).copy(tractions(*it, ghost_type));
    opening_old(*it, ghost_type).copy(opening(*it, ghost_type));
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::assembleResidual(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model->getSpatialDimension();

  Vector<Real> & residual = const_cast<Vector<Real> &>(model->getResidual());

  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type, _ek_cohesive);
  for(; it != last_type; ++it) {
    const Vector<Real> & shapes = fem_cohesive->getShapes(*it, ghost_type);

    Vector<UInt> & elem_filter = element_filter(*it, ghost_type);
    Vector<Real> & traction = tractions(*it, ghost_type);

    UInt size_of_shapes       = shapes.getNbComponent();
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
    UInt nb_quadrature_points = fem_cohesive->getNbQuadraturePoints(*it, ghost_type);

    UInt nb_element = elem_filter.getSize();

    /// compute @f$t_i N_a@f$

    Real * shapes_val       = shapes.storage();
    UInt * elem_filter_val  = elem_filter.storage();

    Vector<Real> * shapes_filtered =
      new Vector<Real>(nb_element*nb_quadrature_points, size_of_shapes, "filtered shapes");
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

    Vector<Real> * traction_cpy = new Vector<Real>(traction);
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
    Vector<Real> * int_t_N = new Vector<Real>(nb_element, spatial_dimension*size_of_shapes,
					      "int_t_N");

    fem_cohesive->integrate(*traction_cpy, *int_t_N,
			    spatial_dimension*size_of_shapes,
			    *it, ghost_type,
			    &elem_filter);

    delete traction_cpy;

    int_t_N->extendComponentsInterlaced(2, int_t_N->getNbComponent() );

    Real * int_t_N_val = int_t_N->storage();
    for (UInt el = 0; el < nb_element; ++el) {
      for (UInt n = 0; n < size_of_shapes*spatial_dimension; ++n)
	int_t_N_val[n] *= -1.;
      int_t_N_val += nb_nodes_per_element*spatial_dimension;
    }

    /// assemble
    model->getFEMBoundary().assembleVector(*int_t_N, residual,
					   model->getDOFSynchronizer().getLocalDOFEquationNumbers(),
					   residual.getNbComponent(),
					   *it, ghost_type, &elem_filter, 1);

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

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type, _ek_cohesive);

  for(; it != last_type; ++it) {
    UInt nb_quadrature_points = fem_cohesive->getNbQuadraturePoints(*it, ghost_type);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);

    const Vector<Real> & shapes = fem_cohesive->getShapes(*it, ghost_type);
    Vector<UInt> & elem_filter = element_filter(*it, ghost_type);

    UInt nb_element = elem_filter.getSize();
    UInt size_of_shapes       = shapes.getNbComponent();

    UInt * elem_filter_it  = elem_filter.storage();

    Vector<Real> * shapes_filtered =
      new Vector<Real>(nb_element*nb_quadrature_points, size_of_shapes, "filtered shapes");

    Vector<Real>::iterator<types::RMatrix> shapes_filtered_it =
      shapes_filtered->begin_reinterpret(size_of_shapes, nb_quadrature_points,
					 nb_element);

    Vector<Real>::const_iterator<types::RMatrix> shapes_it =
      shapes.begin_reinterpret(size_of_shapes, nb_quadrature_points,
			       mesh.getNbElement(*it, ghost_type));

    for (UInt el = 0; el < nb_element; ++el) {
      *shapes_filtered_it = shapes_it[elem_filter_it[el]];
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
    types::RMatrix A(spatial_dimension*size_of_shapes, spatial_dimension*nb_nodes_per_element);

    for ( UInt i = 0; i < spatial_dimension*size_of_shapes; ++i) {
      A(i, i);
      A(i, i + spatial_dimension*size_of_shapes) = -1;
    }

    /// compute traction
    computeTraction(ghost_type);

    /// get the tangent matrix @f$\frac{\partial{(t/\delta)}}{\partial{\delta}} @f$
    Vector<Real> * tangent_stiffness_matrix =
      new Vector<Real>(nb_element*nb_quadrature_points, spatial_dimension*
		       spatial_dimension, "tangent_stiffness_matrix");

    Vector<Real> * normal = new Vector<Real>(nb_element * nb_quadrature_points,
					     spatial_dimension,
					     "normal");
    computeNormal(model->getDisplacement(), *normal, *it, ghost_type);

    tangent_stiffness_matrix->clear();

    computeTangentTraction(*it, *tangent_stiffness_matrix, *normal, ghost_type);

    delete normal;

    UInt size_at_nt_d_n_a = spatial_dimension*nb_nodes_per_element*spatial_dimension*nb_nodes_per_element;
    Vector<Real> * at_nt_d_n_a = new Vector<Real> (nb_element*nb_quadrature_points,
						   size_at_nt_d_n_a,
						   "A^t*N^t*D*N*A");

    Vector<Real>::iterator<types::Vector<Real> > shapes_filt_it = shapes_filtered->begin(size_of_shapes);
    Vector<Real>::iterator<types::RMatrix> D_it = tangent_stiffness_matrix->begin(spatial_dimension, spatial_dimension);
    Vector<Real>::iterator<types::RMatrix> At_Nt_D_N_A_it  = at_nt_d_n_a->begin(spatial_dimension * nb_nodes_per_element,
									      spatial_dimension * nb_nodes_per_element);
    Vector<Real>::iterator<types::RMatrix> At_Nt_D_N_A_end = at_nt_d_n_a->end  (spatial_dimension * nb_nodes_per_element,
									       spatial_dimension * nb_nodes_per_element);

    types::RMatrix N    (spatial_dimension, spatial_dimension * size_of_shapes);
    types::RMatrix N_A  (spatial_dimension, spatial_dimension * nb_nodes_per_element);
    types::RMatrix D_N_A(spatial_dimension, spatial_dimension * nb_nodes_per_element);

    for(; At_Nt_D_N_A_it != At_Nt_D_N_A_end; ++At_Nt_D_N_A_it, ++D_it, ++shapes_filt_it) {
      types::RMatrix & D = *D_it;
      types::RMatrix & At_Nt_D_N_A = *At_Nt_D_N_A_it;
      types::Vector<Real> & shapes_fil = *shapes_filt_it;

      N.clear();
      /**
       * store  the   shapes  in  voigt   notations  matrix  @f$\mathbf{N}  =
       * \begin{array}{cccccc} N_0(\xi) & 0 & N_1(\xi)  &0 & N_2(\xi) & 0 \\
       * 0 & * N_0(\xi)& 0 &N_1(\xi)& 0 & N_2(\xi) \end{array} @f$
       **/
      for (UInt i = 0; i < spatial_dimension ; ++i)
	for (UInt n = 0; n < size_of_shapes; ++n)
	  N(i, i + n) = shapes_fil(n);

      /**
       * compute stiffness matrix  @f$   \mathbf{K}    =    \delta    \mathbf{U}^T
       * \int_{\Gamma_c}    {\mathbf{P}^t \frac{\partial{\mathbf{t}}} {\partial{\delta}}
       * \mathbf{P} d\Gamma \Delta \mathbf{U}}  @f$
       **/
      N_A.mul<false, false>(N, A);
      D_N_A.mul<false, false>(D, N_A);
      At_Nt_D_N_A.mul<true, false>(N_A, D_N_A);
    }

    delete tangent_stiffness_matrix;
    delete shapes_filtered;

    Vector<Real> * K_e = new Vector<Real>(nb_element, size_at_nt_d_n_a,
					  "K_e");

    fem_cohesive->integrate(*at_nt_d_n_a, *K_e,
			    size_at_nt_d_n_a,
			    *it, ghost_type,
			    &elem_filter);

    delete at_nt_d_n_a;

    model->getFEM().assembleMatrix(*K_e, K, spatial_dimension, *it, ghost_type, &elem_filter);
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
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type, _ek_cohesive);

  for(; it != last_type; ++it) {
    Vector<UInt> & elem_filter = element_filter(*it, ghost_type);

    UInt nb_element = elem_filter.getSize();
    UInt nb_quadrature_points = nb_element*fem_cohesive->getNbQuadraturePoints(*it, ghost_type);

    Vector<Real> normal(nb_quadrature_points, spatial_dimension,
			"normal");

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
void MaterialCohesive::computeNormal(const Vector<Real> & position,
				     Vector<Real> & normal,
				     ElementType type,
				     GhostType ghost_type) {
  AKANTU_DEBUG_IN();

#define COMPUTE_NORMAL(type)						\
  fem_cohesive->getShapeFunctions().					\
    computeNormalsOnControlPoints<type, CohesiveReduceFunctionMean>(position, \
								    normal, \
								    ghost_type,	\
								    &(element_filter(type, ghost_type)))

  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_NORMAL);
#undef COMPUTE_NORMAL

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void MaterialCohesive::computeOpening(const Vector<Real> & displacement,
				      Vector<Real> & opening,
				      ElementType type,
				      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

#define COMPUTE_OPENING(type)						\
  fem_cohesive->getShapeFunctions().					\
    interpolateOnControlPoints<type, CohesiveReduceFunctionOpening>(displacement, \
								    opening, \
								    spatial_dimension, \
								    ghost_type, \
								    &(element_filter(type, ghost_type)))

  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_OPENING);
#undef COMPUTE_OPENING

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::computeEnergies() {
  AKANTU_DEBUG_IN();

  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _not_ghost, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _not_ghost, _ek_cohesive);

  Real * memory_space = new Real[2*spatial_dimension];
  types::RVector b(memory_space, spatial_dimension);
  types::RVector h(memory_space + spatial_dimension, spatial_dimension);

  for(; it != last_type; ++it) {
    Vector<Real>::iterator<Real> erev =
      reversible_energy(*it, _not_ghost).begin();
    Vector<Real>::iterator<Real> etot =
      total_energy(*it, _not_ghost).begin();
    Vector<Real>::iterator<types::RVector> traction_it =
      tractions(*it, _not_ghost).begin(spatial_dimension);
    Vector<Real>::iterator<types::RVector> traction_old_it =
      tractions_old(*it, _not_ghost).begin(spatial_dimension);
    Vector<Real>::iterator<types::RVector> opening_it =
      opening(*it, _not_ghost).begin(spatial_dimension);
    Vector<Real>::iterator<types::RVector> opening_old_it =
      opening_old(*it, _not_ghost).begin(spatial_dimension);

    Vector<Real>::iterator<types::RVector>traction_end =
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

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getReversibleEnergy() {
  AKANTU_DEBUG_IN();
  Real erev = 0.;

  /// integrate the dissipated energy for each type of elements
  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _not_ghost, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _not_ghost, _ek_cohesive);

  for(; it != last_type; ++it) {
    erev += fem_cohesive->integrate(reversible_energy(*it, _not_ghost), *it,
				    _not_ghost, &element_filter(*it, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return erev;
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getDissipatedEnergy() {
  AKANTU_DEBUG_IN();
  Real edis = 0.;

  /// integrate the dissipated energy for each type of elements
  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _not_ghost, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _not_ghost, _ek_cohesive);

  for(; it != last_type; ++it) {
    Vector<Real> dissipated_energy(total_energy(*it, _not_ghost));
    dissipated_energy -= reversible_energy(*it, _not_ghost);
    edis += fem_cohesive->integrate(dissipated_energy, *it,
				    _not_ghost, &element_filter(*it, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return edis;
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getEnergy(std::string type) {
  AKANTU_DEBUG_IN();

  if (type == "reversible") return getReversibleEnergy();
  else if (type == "dissipated") return getDissipatedEnergy();

  AKANTU_DEBUG_OUT();
  return 0.;
}


__END_AKANTU__
