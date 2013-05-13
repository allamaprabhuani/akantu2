/**
 * @file   material.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Jul 27 18:15:37 2010
 *
 * @brief  Implementation of the common part of the material class
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
#include "material.hh"
#include "solid_mechanics_model.hh"
#include "sparse_matrix.hh"
#include "dof_synchronizer.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Material::Material(SolidMechanicsModel & model, const ID & id) :
  Memory(model.getMemoryID()),
  id(id),
  name(""),
  model(&model),
  stress("stress", id),
  strain("strain", id),
  element_filter("element_filter", id),
  //  potential_energy_vector(false),
  potential_energy("potential_energy", id),
  is_non_local(false),
  interpolation_inverse_coordinates("interpolation inverse coordinates", id),
  interpolation_points_matrices("interpolation points matrices", id),
  is_init(false) {
  AKANTU_DEBUG_IN();

  finite_deformation=false;

  registerParam("rho",  rho,                 0., ParamAccessType(_pat_parsable | _pat_modifiable), "Density");
  registerParam("id",   this->id,                _pat_readable);
  registerParam("name", name,     std::string(), ParamAccessType(_pat_parsable | _pat_readable));

  spatial_dimension = this->model->getSpatialDimension();

  /// allocate strain stress for local elements
  initInternalArray(strain, spatial_dimension * spatial_dimension);
  initInternalArray(stress, spatial_dimension * spatial_dimension);

  /// for each connectivity types allocate the element filer array of the material
  initInternalArray(element_filter, 1, true);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Material::~Material() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool Material::parseParam(const std::string & key, const std::string & value,
			  __attribute__ ((unused)) const ID & id) {
  try {
    params.parseParam(key, value);
  } catch(...) { return false; }
  return true;
}
/* -------------------------------------------------------------------------- */
void Material::initMaterial() {
  AKANTU_DEBUG_IN();

  resizeInternalArray(stress);
  resizeInternalArray(strain);

  is_init = true;

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template<typename T>
void Material::initInternalArray(ByElementTypeArray<T> & vect,
				 UInt nb_component,
				 bool temporary,
				 ElementKind element_kind) {
  AKANTU_DEBUG_IN();

  model->getFEM().getMesh().initByElementTypeArray(vect,
                                                   nb_component,
                                                   spatial_dimension,
                                                   false,
                                                   element_kind);

  if(!temporary)
    registerInternal(vect);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<> void Material::registerInternal<Real>(ByElementTypeArray<Real> & vect) {
  internal_vectors_real[vect.getID()] = &vect;
}
template<> void Material::registerInternal<UInt>(ByElementTypeArray<UInt> & vect) {
  internal_vectors_uint[vect.getID()] = &vect;
}


/* -------------------------------------------------------------------------- */
template<typename T>
void Material::resizeInternalArray(ByElementTypeArray<T> & by_el_type_vect,
				   ElementKind element_kind) const {
  AKANTU_DEBUG_IN();

  FEM * fem = & model->getFEM();

  if (element_kind == _ek_cohesive)
    fem = & model->getFEM("CohesiveFEM");

  const Mesh & mesh = fem->getMesh();
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    Mesh::type_iterator it  = mesh.firstType(spatial_dimension, gt, element_kind);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, element_kind);
    for(; it != end; ++it) {
      const Array<UInt> & elem_filter = element_filter(*it, gt);

      UInt nb_element           = elem_filter.getSize();
      UInt nb_quadrature_points = fem->getNbQuadraturePoints(*it, gt);
      UInt new_size = nb_element * nb_quadrature_points;

      Array<T> & vect = by_el_type_vect(*it, gt);
      vect.resize(new_size);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute  the  residual  by  assembling  @f$\int_{e}  \sigma_e  \frac{\partial
 * \varphi}{\partial X} dX @f$
 *
 * @param[in] displacements nodes displacements
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void Material::updateResidual(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  computeAllStresses(ghost_type);

  assembleResidual(ghost_type);

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void Material::assembleResidual(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model->getSpatialDimension();

  Array<Real> & residual = const_cast<Array<Real> &>(model->getResidual());

  Mesh & mesh = model->getFEM().getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    const Array<Real> & shapes_derivatives = model->getFEM().getShapesDerivatives(*it, ghost_type);
    Array<UInt> & elem_filter = element_filter(*it, ghost_type);

    UInt size_of_shapes_derivatives = shapes_derivatives.getNbComponent();
    UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(*it);
    UInt nb_quadrature_points       = model->getFEM().getNbQuadraturePoints(*it, ghost_type);

    UInt nb_element = elem_filter.getSize();

    /// compute @f$\sigma \frac{\partial \varphi}{\partial X}@f$ by @f$\mathbf{B}^t \mathbf{\sigma}_q@f$
    Array<Real> * sigma_dphi_dx =
      new Array<Real>(nb_element*nb_quadrature_points,
		      size_of_shapes_derivatives, "sigma_x_dphi_/_dX");

    Array<Real> * shapesd_filtered =
      new Array<Real>(0, size_of_shapes_derivatives, "filtered shapesd");

    FEM::filterElementalData(mesh, shapes_derivatives, *shapesd_filtered,
                             *it, ghost_type, elem_filter);

    Array<Real> & stress_vect = stress(*it, ghost_type);

    Array<Real>::iterator< Matrix<Real> > sigma =
      stress_vect.begin(spatial_dimension, spatial_dimension);
    Array<Real>::iterator< Matrix<Real> > B =
      shapesd_filtered->begin(spatial_dimension, nb_nodes_per_element);
    Array<Real>::iterator< Matrix<Real> > Bt_sigma_it =
      sigma_dphi_dx->begin(spatial_dimension, nb_nodes_per_element);

    for (UInt q = 0; q < nb_element*nb_quadrature_points; ++q, ++sigma, ++B, ++Bt_sigma_it)
      Bt_sigma_it->mul<false,false>(*sigma, *B);

    delete shapesd_filtered;

    /**
     * compute @f$\int \sigma  * \frac{\partial \varphi}{\partial X}dX@f$ by  @f$ \sum_q \mathbf{B}^t
     * \mathbf{\sigma}_q \overline w_q J_q@f$
     */
    Array<Real> * int_sigma_dphi_dx = new Array<Real>(nb_element,
						      nb_nodes_per_element * spatial_dimension,
						      "int_sigma_x_dphi_/_dX");

    model->getFEM().integrate(*sigma_dphi_dx, *int_sigma_dphi_dx,
                              size_of_shapes_derivatives,
                              *it, ghost_type,
                              elem_filter);
    delete sigma_dphi_dx;

    /// assemble
    model->getFEM().assembleArray(*int_sigma_dphi_dx, residual,
                                  model->getDOFSynchronizer().getLocalDOFEquationNumbers(),
                                  residual.getNbComponent(),
                                  *it, ghost_type, elem_filter, -1);
    delete int_sigma_dphi_dx;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute  the  stress from the strain
 *
 * @param[in] current_position nodes postition + displacements
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void Material::computeAllStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  resizeInternalArray(stress);
  resizeInternalArray(strain);

  UInt spatial_dimension = model->getSpatialDimension();

  Mesh::type_iterator it = model->getFEM().getMesh().firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = model->getFEM().getMesh().lastType(spatial_dimension, ghost_type);

  for(; it != last_type; ++it) {
    Array<UInt> & elem_filter = element_filter(*it, ghost_type);
    Array<Real> & strain_vect = strain(*it, ghost_type);

    /// compute @f$\nabla u@f$
    model->getFEM().gradientOnQuadraturePoints(model->getDisplacement(), strain_vect,
                                               spatial_dimension,
                                               *it, ghost_type, elem_filter);

    /// compute @f$\mathbf{\sigma}_q@f$ from @f$\nabla u@f$
    computeStress(*it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::setToSteadyState(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const Array<Real> & displacement = model->getDisplacement();

  resizeInternalArray(strain);

  UInt spatial_dimension = model->getSpatialDimension();

  Mesh::type_iterator it = model->getFEM().getMesh().firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = model->getFEM().getMesh().lastType(spatial_dimension, ghost_type);

  for(; it != last_type; ++it) {
    Array<UInt> & elem_filter = element_filter(*it, ghost_type);
    Array<Real> & strain_vect = strain(*it, ghost_type);

    /// compute @f$\nabla u@f$
    model->getFEM().gradientOnQuadraturePoints(displacement, strain_vect,
					       spatial_dimension,
					       *it, ghost_type, elem_filter);

    setToSteadyState(*it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/**
 * Compute  the stiffness  matrix by  assembling @f$\int_{\omega}  B^t  \times D
 * \times B d\omega @f$
 *
 * @param[in] current_position nodes postition + displacements
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void Material::assembleStiffnessMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model->getSpatialDimension();

  Mesh & mesh = model->getFEM().getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    switch(spatial_dimension) {
    case 1: { assembleStiffnessMatrix<1>(*it, ghost_type); break; }
    case 2: { assembleStiffnessMatrix<2>(*it, ghost_type); break; }
    case 3: { assembleStiffnessMatrix<3>(*it, ghost_type); break; }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void Material::assembleStiffnessMatrix(const ElementType & type,
                                       GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  SparseMatrix & K = const_cast<SparseMatrix &>(model->getStiffnessMatrix());

  const Array<Real> & shapes_derivatives = model->getFEM().getShapesDerivatives(type,ghost_type);

  Array<UInt> & elem_filter = element_filter(type, ghost_type);
  Array<Real> & strain_vect = strain(type, ghost_type);

  UInt nb_element                 = elem_filter.getSize();
  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points       = model->getFEM().getNbQuadraturePoints(type, ghost_type);

  strain_vect.resize(nb_quadrature_points * nb_element);

  model->getFEM().gradientOnQuadraturePoints(model->getDisplacement(), strain_vect,
                                             dim, type, ghost_type, elem_filter);

  UInt tangent_size = getTangentStiffnessVoigtSize(dim);

  Array<Real> * tangent_stiffness_matrix =
    new Array<Real>(nb_element*nb_quadrature_points, tangent_size * tangent_size,
		    "tangent_stiffness_matrix");

  tangent_stiffness_matrix->clear();

  computeTangentModuli(type, *tangent_stiffness_matrix, ghost_type);

  Array<Real> * shapesd_filtered =
    new Array<Real>(0, dim * nb_nodes_per_element, "filtered shapesd");

  FEM::filterElementalData(model->getFEM().getMesh(), shapes_derivatives, *shapesd_filtered,
                           type, ghost_type, elem_filter);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  UInt bt_d_b_size = dim * nb_nodes_per_element;

  Array<Real> * bt_d_b = new Array<Real>(nb_element * nb_quadrature_points,
					 bt_d_b_size * bt_d_b_size,
					 "B^t*D*B");

  Matrix<Real> B(tangent_size, dim * nb_nodes_per_element);
  Matrix<Real> Bt_D(dim * nb_nodes_per_element, tangent_size);

  Array<Real>::iterator< Matrix<Real> > shapes_derivatives_filtered_it =
    shapesd_filtered->begin(dim, nb_nodes_per_element);

  Array<Real>::iterator< Matrix<Real> > Bt_D_B_it = bt_d_b->begin(dim*nb_nodes_per_element,
                                                                  dim*nb_nodes_per_element);

  Array<Real>::iterator< Matrix<Real> > D_it  = tangent_stiffness_matrix->begin(tangent_size,
                                                                                tangent_size);
  Array<Real>::iterator< Matrix<Real> > D_end = tangent_stiffness_matrix->end  (tangent_size,
                                                                                tangent_size);


  for(; D_it != D_end; ++D_it, ++Bt_D_B_it, ++shapes_derivatives_filtered_it) {
    Matrix<Real> & D = *D_it;
    Matrix<Real> & Bt_D_B = *Bt_D_B_it;

    transferBMatrixToSymVoigtBMatrix<dim>(*shapes_derivatives_filtered_it, B, nb_nodes_per_element);
    Bt_D.mul<true, false>(B, D);
    Bt_D_B.mul<false, false>(Bt_D, B);
  }

  delete tangent_stiffness_matrix;
  delete shapesd_filtered;

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  Array<Real> * K_e = new Array<Real>(nb_element,
				      bt_d_b_size * bt_d_b_size,
				      "K_e");

  model->getFEM().integrate(*bt_d_b, *K_e,
                            bt_d_b_size * bt_d_b_size,
                            type, ghost_type,
                            elem_filter);

  delete bt_d_b;

  model->getFEM().assembleMatrix(*K_e, K, spatial_dimension, type, ghost_type, elem_filter);
  delete K_e;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::computeAllStressesFromTangentModuli(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model->getSpatialDimension();

  Mesh & mesh = model->getFEM().getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    switch(spatial_dimension) {
    case 1: { computeAllStressesFromTangentModuli<1>(*it, ghost_type); break; }
    case 2: { computeAllStressesFromTangentModuli<2>(*it, ghost_type); break; }
    case 3: { computeAllStressesFromTangentModuli<3>(*it, ghost_type); break; }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void Material::computeAllStressesFromTangentModuli(const ElementType & type,
                                                   GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const Array<Real> & shapes_derivatives = model->getFEM().getShapesDerivatives(type, ghost_type);
  Array<UInt> & elem_filter = element_filter(type, ghost_type);
  Array<Real> & strain_vect = strain(type, ghost_type);

  UInt nb_element                 = elem_filter.getSize();
  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points       = model->getFEM().getNbQuadraturePoints(type, ghost_type);


  strain_vect.resize(nb_quadrature_points * nb_element);

  Array<Real> & disp = model->getDisplacement();

  model->getFEM().gradientOnQuadraturePoints(disp, strain_vect,
                                             dim, type, ghost_type, elem_filter);

  UInt tangent_moduli_size = getTangentStiffnessVoigtSize(dim);

  Array<Real> * tangent_moduli_tensors =
    new Array<Real>(nb_element*nb_quadrature_points, tangent_moduli_size * tangent_moduli_size,
		    "tangent_moduli_tensors");

  tangent_moduli_tensors->clear();
  computeTangentModuli(type, *tangent_moduli_tensors, ghost_type);

  Array<Real> * shapesd_filtered =
    new Array<Real>(0, dim* nb_nodes_per_element, "filtered shapesd");

  FEM::filterElementalData(model->getFEM().getMesh(), shapes_derivatives, *shapesd_filtered,
			   type, ghost_type, elem_filter);

  Array<Real> filtered_u(nb_element, nb_nodes_per_element * spatial_dimension);

  FEM::extractNodalToElementField(model->getFEM().getMesh(), disp, filtered_u,
                                  type, ghost_type, elem_filter);

  /// compute @f$\mathbf{D} \mathbf{B} \mathbf{u}@f$
  Array<Real>::iterator< Matrix<Real> > shapes_derivatives_filtered_it =
    shapesd_filtered->begin(dim, nb_nodes_per_element);

  Array<Real>::iterator< Matrix<Real> > D_it  = tangent_moduli_tensors->begin(tangent_moduli_size,
									      tangent_moduli_size);
  Array<Real>::iterator< Matrix<Real> > sigma_it  = stress(type, ghost_type).begin(spatial_dimension,
										   spatial_dimension);
  Array<Real>::iterator< Vector<Real> > u_it = filtered_u.begin(spatial_dimension * nb_nodes_per_element);

  Matrix<Real> B(tangent_moduli_size, spatial_dimension * nb_nodes_per_element);
  Vector<Real> Bu(tangent_moduli_size);
  Vector<Real> DBu(tangent_moduli_size);

  for (UInt e = 0; e < nb_element; ++e, ++u_it) {
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++D_it, ++shapes_derivatives_filtered_it, ++sigma_it) {
      Vector<Real> & u = *u_it;
      Matrix<Real> & sigma = *sigma_it;
      Matrix<Real> & D = *D_it;

      transferBMatrixToSymVoigtBMatrix<dim>(*shapes_derivatives_filtered_it, B, nb_nodes_per_element);

      Bu.mul<false>(B, u);
      DBu.mul<false>(D, Bu);

      // Voigt notation to full symmetric tensor
      for (UInt i = 0; i < dim; ++i) sigma(i, i) = DBu(i);
      if(dim == 2) {
        sigma(0,1) = sigma(1,0) = DBu(2);
      } else if(dim == 3) {
        sigma(1,2) = sigma(2,1) = DBu(3);
        sigma(0,2) = sigma(2,0) = DBu(4);
        sigma(0,1) = sigma(1,0) = DBu(5);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}



/* -------------------------------------------------------------------------- */
void Material::computePotentialEnergy(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if(ghost_type != _not_ghost) return;
  Real * epot = potential_energy(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  computePotentialEnergyOnQuad(grad_u, sigma, *epot);
  epot++;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::computePotentialEnergyByElement() {
  AKANTU_DEBUG_IN();

  Mesh & mesh = model->getFEM().getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension);
  for(; it != last_type; ++it) {

    if(!potential_energy.exists(*it, _not_ghost)) {
      UInt nb_element = element_filter(*it, _not_ghost).getSize();
      UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(*it, _not_ghost);

      potential_energy.alloc(nb_element * nb_quadrature_points, 1,
                             *it, _not_ghost);
    }

    computePotentialEnergy(*it);
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */

void Material::computePotentialEnergyByElement(ElementType type, UInt index,
                                               Vector<Real> & epot_on_quad_points){

  Array<Real>::iterator< Matrix<Real> > strain_it =
    this->strain(type).begin(spatial_dimension,
                             spatial_dimension);
  Array<Real>::iterator< Matrix<Real> > strain_end =
    this->strain(type).begin(spatial_dimension,
                             spatial_dimension);
  Array<Real>::iterator< Matrix<Real> > stress_it =
    this->stress(type).begin(spatial_dimension,
                             spatial_dimension);

  UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(type);

  strain_it  += index*nb_quadrature_points;
  strain_end += (index+1)*nb_quadrature_points;
  stress_it  += index*nb_quadrature_points;

  Real * epot_quad = epot_on_quad_points.storage();

  for(;strain_it != strain_end; ++strain_it, ++stress_it, ++epot_quad) {
    Matrix<Real> & grad_u = *strain_it;
    Matrix<Real> & sigma  = *stress_it;
    computePotentialEnergyOnQuad(grad_u, sigma, *epot_quad);
  }
}

/* -------------------------------------------------------------------------- */
Real Material::getPotentialEnergy() {
  AKANTU_DEBUG_IN();
  Real epot = 0.;

  computePotentialEnergyByElement();

  /// integrate the potential energy for each type of elements
  Mesh & mesh = model->getFEM().getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension);
  for(; it != last_type; ++it) {

    epot += model->getFEM().integrate(potential_energy(*it, _not_ghost), *it,
                                      _not_ghost, element_filter(*it, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */
Real Material::getPotentialEnergy(ElementType & type, UInt index) {
  AKANTU_DEBUG_IN();
  Real epot = 0.;

  Vector<Real> epot_on_quad_points(model->getFEM().getNbQuadraturePoints(type));

  computePotentialEnergyByElement(type,index,epot_on_quad_points);

  epot = model->getFEM().integrate(epot_on_quad_points, type, element_filter(type)(index));

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */
Real Material::getEnergy(std::string type) {
  AKANTU_DEBUG_IN();
  if(type == "potential") return getPotentialEnergy();
  AKANTU_DEBUG_OUT();
  return 0.;
}

/* -------------------------------------------------------------------------- */
Real Material::getEnergy(std::string energy_id, ElementType type, UInt index) {
  AKANTU_DEBUG_IN();
  if(energy_id == "potential") return getPotentialEnergy(type, index);
  AKANTU_DEBUG_OUT();
  return 0.;
}

/* -------------------------------------------------------------------------- */
void Material::computeQuadraturePointsCoordinates(ByElementTypeReal & quadrature_points_coordinates, const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  const Mesh & mesh = model->getFEM().getMesh();

  Array<Real> nodes_coordinates(mesh.getNodes(), true);
  nodes_coordinates += model->getDisplacement();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    const Array<UInt> & elem_filter = element_filter(*it, ghost_type);

    UInt nb_element  = elem_filter.getSize();
    UInt nb_tot_quad = model->getFEM().getNbQuadraturePoints(*it, ghost_type) * nb_element;

    Array<Real> & quads = quadrature_points_coordinates(*it, ghost_type);
    quads.resize(nb_tot_quad);

    model->getFEM().interpolateOnQuadraturePoints(nodes_coordinates,
                                                  quads, spatial_dimension,
                                                  *it, ghost_type, elem_filter);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::initElementalFieldInterpolation(ByElementTypeReal & interpolation_points_coordinates) {
  AKANTU_DEBUG_IN();
  const Mesh & mesh = model->getFEM().getMesh();

  ByElementTypeReal quadrature_points_coordinates("quadrature_points_coordinates_tmp_nl", id);
  mesh.initByElementTypeArray(quadrature_points_coordinates,
			      spatial_dimension, spatial_dimension);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType ghost_type = *gt;

    computeQuadraturePointsCoordinates(quadrature_points_coordinates, ghost_type);

    Mesh::type_iterator it   = mesh.firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last = mesh.lastType(spatial_dimension, ghost_type);
    for (; it != last; ++it) {
      ElementType type = *it;
      UInt nb_element = mesh.getNbElement(type, ghost_type);
      if (nb_element == 0) continue;

#define AKANTU_INIT_INTERPOLATE_ELEMENTAL_FIELD(type)			\
      initElementalFieldInterpolation<type>(quadrature_points_coordinates(type, ghost_type), \
                                            interpolation_points_coordinates(type, ghost_type),	\
                                            ghost_type)			\

      AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(AKANTU_INIT_INTERPOLATE_ELEMENTAL_FIELD);
#undef AKANTU_INIT_INTERPOLATE_ELEMENTAL_FIELD
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void Material::initElementalFieldInterpolation(const Array<Real> & quad_coordinates,
                                               const Array<Real> & interpolation_points_coordinates,
                                               const GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  UInt size_inverse_coords =
    getSizeElementalFieldInterpolationCoodinates<type>(ghost_type);

  Array<UInt> & elem_fil = element_filter(type, ghost_type);
  UInt nb_element = elem_fil.getSize();
  UInt nb_quad_per_element = model->getFEM().getNbQuadraturePoints(type, ghost_type);

  UInt nb_interpolation_points = interpolation_points_coordinates.getSize();
  AKANTU_DEBUG_ASSERT(nb_interpolation_points % nb_element == 0,
                      "Can't interpolate elemental field on elements, the coordinates vector has a wrong size");
  UInt nb_interpolation_points_per_elem = nb_interpolation_points / nb_element;

  if(!interpolation_inverse_coordinates.exists(type, ghost_type))
    interpolation_inverse_coordinates.alloc(nb_element,
                                            size_inverse_coords*size_inverse_coords,
                                            type, ghost_type);

  if(!interpolation_points_matrices.exists(type, ghost_type))
    interpolation_points_matrices.alloc(nb_element,
                                        nb_interpolation_points_per_elem * size_inverse_coords,
                                        type, ghost_type);

  Array<Real> & interp_inv_coord = interpolation_inverse_coordinates(type, ghost_type);
  Array<Real> & interp_points_mat = interpolation_points_matrices(type, ghost_type);

  Matrix<Real> quad_coord_matrix(size_inverse_coords, size_inverse_coords);

  Array<Real>::const_iterator< Matrix<Real> > quad_coords_it =
    quad_coordinates.begin_reinterpret(spatial_dimension,
                                       nb_quad_per_element,
                                       nb_element);

  Array<Real>::const_iterator< Matrix<Real> > points_coords_it =
    interpolation_points_coordinates.begin_reinterpret(spatial_dimension,
                                                       nb_interpolation_points_per_elem,
                                                       nb_element);

  Array<Real>::iterator< Matrix<Real> > inv_quad_coord_it =
    interp_inv_coord.begin(size_inverse_coords, size_inverse_coords);

  Array<Real>::iterator< Matrix<Real> > inv_points_mat_it =
    interp_points_mat.begin(nb_interpolation_points_per_elem, size_inverse_coords);

  /// loop over the elements of the current material and element type
  for (UInt el = 0; el < nb_element; ++el) {
    /// matrix containing the quadrature points coordinates
    const Matrix<Real> & quad_coords = *quad_coords_it;
    /// matrix to store the matrix inversion result
    Matrix<Real> & inv_quad_coord_matrix = *inv_quad_coord_it;

    /// insert the quad coordinates in a matrix compatible with the interpolation
    buildElementalFieldInterpolationCoodinates<type>(quad_coords,
                                                     quad_coord_matrix);

    /// invert the interpolation matrix
    inv_quad_coord_matrix.inverse(quad_coord_matrix);


    /// matrix containing the interpolation points coordinates
    const Matrix<Real> & points_coords = *points_coords_it;
    /// matrix to store the interpolation points coordinates
    /// compatible with these functions
    Matrix<Real> & inv_points_coord_matrix = *inv_points_mat_it;

    /// insert the quad coordinates in a matrix compatible with the interpolation
    buildElementalFieldInterpolationCoodinates<type>(points_coords,
                                                     inv_points_coord_matrix);


    ++inv_quad_coord_it;
    ++inv_points_mat_it;
    ++quad_coords_it;
    ++points_coords_it;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::interpolateStress(ByElementTypeReal & result) {
  AKANTU_DEBUG_IN();

  const Mesh & mesh = model->getFEM().getMesh();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType ghost_type = *gt;

    Mesh::type_iterator it   = mesh.firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last = mesh.lastType(spatial_dimension, ghost_type);
    for (; it != last; ++it) {
      ElementType type = *it;
      UInt nb_element = mesh.getNbElement(type, ghost_type);
      if (nb_element == 0) continue;

      Array<Real> & res = result(type, ghost_type);

#define INTERPOLATE_ELEMENTAL_FIELD(type)			\
      interpolateElementalField<type>(stress(type, ghost_type),	\
                                      res,			\
                                      ghost_type)

      AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(INTERPOLATE_ELEMENTAL_FIELD);
#undef INTERPOLATE_ELEMENTAL_FIELD
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void Material::interpolateElementalField(const Array<Real> & field,
                                         Array<Real> & result,
                                         const GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<UInt> & elem_fil = element_filter(type, ghost_type);
  UInt nb_element = elem_fil.getSize();
  UInt nb_quad_per_element = model->getFEM().getNbQuadraturePoints(type, ghost_type);
  UInt size_inverse_coords = getSizeElementalFieldInterpolationCoodinates<type>(ghost_type);

  Matrix<Real> coefficients(nb_quad_per_element, field.getNbComponent());

  const Array<Real> & interp_inv_coord = interpolation_inverse_coordinates(type,
                                                                           ghost_type);
  const Array<Real> & interp_points_coord = interpolation_points_matrices(type,
                                                                          ghost_type);

  UInt nb_interpolation_points_per_elem = interp_points_coord.getNbComponent() / size_inverse_coords;

  Array<Real>::const_iterator< Matrix<Real> > field_it
    = field.begin_reinterpret(field.getNbComponent(),
                              nb_quad_per_element,
                              nb_element);

  Array<Real>::const_iterator< Matrix<Real> > interpolation_points_coordinates_it =
    interp_points_coord.begin(nb_interpolation_points_per_elem, size_inverse_coords);

  Array<Real>::iterator< Matrix<Real> > result_it
    = result.begin_reinterpret(field.getNbComponent(),
                               nb_interpolation_points_per_elem,
                               nb_element);

  Array<Real>::const_iterator< Matrix<Real> > inv_quad_coord_it =
    interp_inv_coord.begin(size_inverse_coords, size_inverse_coords);

  /// loop over the elements of the current material and element type
  for (UInt el = 0; el < nb_element;
       ++el, ++field_it, ++result_it,
         ++inv_quad_coord_it, ++interpolation_points_coordinates_it) {
    /**
     * matrix containing the inversion of the quadrature points'
     * coordinates
     */
    const Matrix<Real> & inv_quad_coord_matrix = *inv_quad_coord_it;

    /**
     * multiply it by the field values over quadrature points to get
     * the interpolation coefficients
     */
    coefficients.mul<false, true>(inv_quad_coord_matrix, *field_it);

    /// matrix containing the points' coordinates
    const Matrix<Real> & coord = *interpolation_points_coordinates_it;

    /// multiply the coordinates matrix by the coefficients matrix and store the result
    (*result_it).mul<true, true>(coefficients, coord);
  }


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

const Array<Real> & Material::getArray(const ID & vect_id, const ElementType & type, const GhostType & ghost_type) const {
  std::stringstream sstr;
  std::string ghost_id = "";
  if (ghost_type == _ghost) ghost_id = ":ghost";
  sstr << id << ":" << vect_id << ":" << type << ghost_id;

  ID fvect_id = sstr.str();
  try {
    return Memory::getArray<Real>(fvect_id);
  } catch(debug::Exception & e) {
    AKANTU_EXCEPTION("The material " << name << "(" <<id << ") does not contain a vector " << vect_id << "(" << fvect_id << ") [" << e << "]");
  }
}

/* -------------------------------------------------------------------------- */
Array<Real> & Material::getArray(const ID & vect_id, const ElementType & type, const GhostType & ghost_type) {
  std::stringstream sstr;
  std::string ghost_id = "";
  if (ghost_type == _ghost) ghost_id = ":ghost";
  sstr << id << ":" << vect_id << ":" << type << ghost_id;

  ID fvect_id = sstr.str();
  try {
    return Memory::getArray<Real>(fvect_id);
  } catch(debug::Exception & e) {
    AKANTU_EXCEPTION("The material " << name << "(" <<id << ") does not contain a vector " << vect_id << "(" << fvect_id << ") [" << e << "]");
  }
}

/* -------------------------------------------------------------------------- */
void Material::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  std::string type = id.substr(id.find_last_of(":") + 1);

  stream << space << "Material " << type << " [" << std::endl;
  params.printself(stream, indent);
  stream << space << "]" << std::endl;
}


/* -------------------------------------------------------------------------- */
template void Material::initInternalArray<Real>(ByElementTypeArray<Real> & vect,
						UInt nb_component,
						bool temporary,
						ElementKind element_kind);

template void Material::initInternalArray<UInt>(ByElementTypeArray<UInt> & vect,
						UInt nb_component,
						bool temporary,
						ElementKind element_kind);

template void Material::initInternalArray<Int>(ByElementTypeArray<Int> & vect,
					       UInt nb_component,
					       bool temporary,
					       ElementKind element_kind);

template void Material::initInternalArray<bool>(ByElementTypeArray<bool> & vect,
						UInt nb_component,
						bool temporary,
						ElementKind element_kind);


template void Material::resizeInternalArray<Real>(ByElementTypeArray<Real> & vect,
						  ElementKind element_kind) const;

template void Material::resizeInternalArray<UInt>(ByElementTypeArray<UInt> & vect,
						  ElementKind element_kind) const;

template void Material::resizeInternalArray<Int>(ByElementTypeArray<Int> & vect,
						 ElementKind element_kind) const;

template void Material::resizeInternalArray<bool>(ByElementTypeArray<bool> & vect,
						  ElementKind element_kind) const;


__END_AKANTU__
