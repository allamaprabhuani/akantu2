/**
 * @file   material.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:43:41 2010
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
Material::Material(Model & model, const MaterialID & id) : Memory(model.getMemoryID()),
							   id(id),
							   name(""),
							   potential_energy_vector(false),
							   is_init(false) {
  AKANTU_DEBUG_IN();

  this->model = dynamic_cast<SolidMechanicsModel*>(&model);
  AKANTU_DEBUG_ASSERT(this->model,"model has wrong type: cannot proceed");
  spatial_dimension = this->model->getSpatialDimension();

  /// allocate strain stress for local elements
  initInternalVector(strain, spatial_dimension*spatial_dimension, "strain");
  initInternalVector(stress, spatial_dimension*spatial_dimension, "stress");

  /// for each connectivity types allocate the element filer array of the material
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    std::string ghost_id = "";
    if (gt == _ghost) ghost_id = "ghost_";

    const Mesh::ConnectivityTypeList & type_list =
      model.getFEM().getMesh().getConnectivityTypeList(gt);
    Mesh::ConnectivityTypeList::const_iterator it;
    for(it = type_list.begin(); it != type_list.end(); ++it) {
      if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;
      std::stringstream sstr; sstr << id << ":" << ghost_id << "element_filer:"<< *it;
      element_filter(*it, gt) = &(alloc<UInt> (sstr.str(), 0, 1));
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Material::~Material() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::setParam(const std::string & key, const std::string & value,
			__attribute__ ((unused)) const MaterialID & id) {
  if(key == "name") name = std::string(value);
  else AKANTU_DEBUG_ERROR("Unknown material property : " << key);
}
/* -------------------------------------------------------------------------- */
void Material::initMaterial() {
  AKANTU_DEBUG_IN();
  resizeInternalVector(stress);
  resizeInternalVector(strain);
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template<typename T>
void Material::initInternalVector(ByElementTypeVector<T> & vect,
				  UInt nb_component,
				  const ID & vect_id) {
  AKANTU_DEBUG_IN();

  model->getFEM().getMesh().initByElementTypeVector(vect, nb_component, spatial_dimension,
						    id, vect_id);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template<typename T>
void Material::resizeInternalVector(ByElementTypeVector<T> & by_el_type_vect) {
  AKANTU_DEBUG_IN();

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    const Mesh::ConnectivityTypeList & type_list =
      model->getFEM().getMesh().getConnectivityTypeList(gt);
    Mesh::ConnectivityTypeList::const_iterator it;
    for(it = type_list.begin(); it != type_list.end(); ++it) {
      if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

      Vector<UInt> * elem_filter = element_filter(*it, gt);

      UInt nb_element           = elem_filter->getSize();
      UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(*it, gt);
      UInt new_size = nb_element*nb_quadrature_points;

      Vector<T> * vect = by_el_type_vect(*it, gt);
      UInt size = vect->getSize();
      UInt nb_component = vect->getNbComponent();

      vect->resize(new_size);
      memset(vect->values + size * nb_component, 0, (new_size - size) * nb_component * sizeof(T));
    }
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
/**
 * Compute  the  residual  by  assembling  @f$\int_{e}  \sigma_e  \frac{\partial
 * \varphi}{\partial X} dX @f$
 *
 * @param[in] current_position nodes postition + displacements
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void Material::updateResidual(Vector<Real> & current_position, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model->getSpatialDimension();

  Vector<Real> & residual = const_cast<Vector<Real> &>(model->getResidual());

  const Mesh::ConnectivityTypeList & type_list =
    model->getFEM().getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {

    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    const Vector<Real> & shapes_derivatives = model->getFEM().getShapesDerivatives(*it, ghost_type);
    Vector<UInt> * elem_filter = element_filter(*it, ghost_type);
    Vector<Real> * strain_vect = strain(*it, ghost_type);
    Vector<Real> * stress_vect = stress(*it, ghost_type);


    UInt size_of_shapes_derivatives = shapes_derivatives.getNbComponent();
    UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(*it);
    UInt nb_quadrature_points       = model->getFEM().getNbQuadraturePoints(*it, ghost_type);

    UInt nb_element = elem_filter->getSize();

    /// compute @f$\nabla u@f$
    model->getFEM().gradientOnQuadraturePoints(current_position, *strain_vect,
					      spatial_dimension,
					      *it, ghost_type, elem_filter);

    /// compute @f$\mathbf{\sigma}_q@f$ from @f$\nabla u@f$
    computeStress(*it, ghost_type);

    /// compute @f$\sigma \frac{\partial \varphi}{\partial X}@f$ by @f$\mathbf{B}^t \mathbf{\sigma}_q@f$
    Vector<Real> * sigma_dphi_dx =
      new Vector<Real>(nb_element*nb_quadrature_points, size_of_shapes_derivatives, "sigma_x_dphi_/_dX");

    Real * shapesd           = shapes_derivatives.values;
    UInt size_of_shapesd     = shapes_derivatives.getNbComponent();
    Real * shapesd_val;
    UInt * elem_filter_val   = elem_filter->values;

    Vector<Real> * shapesd_filtered =
      new Vector<Real>(nb_element*nb_quadrature_points, size_of_shapes_derivatives, "filtered shapesd");
    Real * shapesd_filtered_val = shapesd_filtered->values;

    for (UInt el = 0; el < nb_element; ++el) {
      shapesd_val = shapesd + elem_filter_val[el] * size_of_shapesd * nb_quadrature_points;
      memcpy(shapesd_filtered_val,
	     shapesd_val,
	     size_of_shapesd * nb_quadrature_points * sizeof(Real));
      shapesd_filtered_val += size_of_shapesd * nb_quadrature_points;
    }

    Math::matrix_matrixt(nb_nodes_per_element, spatial_dimension, spatial_dimension,
			 *shapesd_filtered,
			 *stress_vect,
			 *sigma_dphi_dx);

    delete shapesd_filtered;

    // for (UInt el = 0; el < nb_element; ++el) {
    //   shapesd_val = shapesd + elem_filter_val[el]*size_of_shapesd*nb_quadrature_points;
    //   for (UInt q = 0; q < nb_quadrature_points; ++q) {
    // 	Math::matrix_matrixt(nb_nodes_per_element, spatial_dimension, spatial_dimension,
    // 			     shapesd_val, stress_val, sigma_dphi_dx_val);
    // 	shapesd_val       += offset_shapesd_val;
    // 	stress_val        += offset_stress_val;
    // 	sigma_dphi_dx_val += offset_sigma_dphi_dx_val;
    //   }
    // }

    /**
     * compute @f$\int \sigma  * \frac{\partial \varphi}{\partial X}dX@f$ by  @f$ \sum_q \mathbf{B}^t
     * \mathbf{\sigma}_q \overline w_q J_q@f$
     */

    Vector<Real> * int_sigma_dphi_dx = new Vector<Real>(0, nb_nodes_per_element * spatial_dimension,
							"int_sigma_x_dphi_/_dX");

    model->getFEM().integrate(*sigma_dphi_dx, *int_sigma_dphi_dx,
			      size_of_shapes_derivatives,
			      *it, ghost_type,
			      elem_filter);
    delete sigma_dphi_dx;

    /// assemble
    model->getFEM().assembleVector(*int_sigma_dphi_dx, residual,
				   model->getDOFSynchronizer().getLocalDOFEquationNumbers(),
				   residual.getNbComponent(),
				   *it, ghost_type, elem_filter, -1);
    delete int_sigma_dphi_dx;
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
void Material::assembleStiffnessMatrix(Vector<Real> & current_position, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model->getSpatialDimension();

  const Mesh::ConnectivityTypeList & type_list =
    model->getFEM().getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {

    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    switch(spatial_dimension) {
    case 1: { assembleStiffnessMatrix<1>(current_position, *it, ghost_type); break; }
    case 2: { assembleStiffnessMatrix<2>(current_position, *it, ghost_type); break; }
    case 3: { assembleStiffnessMatrix<3>(current_position, *it, ghost_type); break; }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void Material::assembleStiffnessMatrix(Vector<Real> & current_position,
				       const ElementType & type,
				       GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  SparseMatrix & K = const_cast<SparseMatrix &>(model->getStiffnessMatrix());

  const Vector<Real> & shapes_derivatives = model->getFEM().getShapesDerivatives(type,ghost_type);

  Vector<UInt> * elem_filter = element_filter(type, ghost_type);
  Vector<Real> * strain_vect = strain(type, ghost_type);

  UInt * elem_filter_val = elem_filter->values;

  UInt nb_element                 = elem_filter->getSize();
  UInt size_of_shapes_derivatives = shapes_derivatives.getNbComponent();
  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points       = model->getFEM().getNbQuadraturePoints(type, ghost_type);

  model->getFEM().gradientOnQuadraturePoints(current_position, *strain_vect,
					     dim, type, ghost_type, elem_filter);

  UInt tangent_size = getTangentStiffnessVoigtSize(dim);

  Vector<Real> * tangent_stiffness_matrix =
    new Vector<Real>(nb_element*nb_quadrature_points, tangent_size * tangent_size,
		     "tangent_stiffness_matrix");

  computeTangentStiffness(type, *tangent_stiffness_matrix, ghost_type);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  UInt bt_d_b_size = dim * nb_nodes_per_element;

  Vector<Real> * bt_d_b = new Vector<Real>(nb_element*nb_quadrature_points,
					   bt_d_b_size * bt_d_b_size,
					   "B^t*D*B");

  UInt size_of_b = tangent_size * bt_d_b_size;
  Real * B = new Real[size_of_b];
  Real * Bt_D = new Real[size_of_b];
  Real * Bt_D_B = bt_d_b->values;
  Real * D = tangent_stiffness_matrix->values;

  UInt offset_bt_d_b = bt_d_b_size * bt_d_b_size;
  UInt offset_d      = tangent_size * tangent_size;

  for (UInt e = 0; e < nb_element; ++e) {
    Real * shapes_derivatives_val =
      shapes_derivatives.values + elem_filter_val[e]*size_of_shapes_derivatives*nb_quadrature_points;

    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      transferBMatrixToSymVoigtBMatrix<dim>(shapes_derivatives_val, B, nb_nodes_per_element);
      Math::matrixt_matrix(bt_d_b_size, tangent_size, tangent_size, B, D, Bt_D);
      Math::matrix_matrix(bt_d_b_size, bt_d_b_size, tangent_size, Bt_D, B, Bt_D_B);

      shapes_derivatives_val += size_of_shapes_derivatives;
      D      += offset_d;
      Bt_D_B += offset_bt_d_b;
    }
  }

  delete [] B;
  delete [] Bt_D;

  delete tangent_stiffness_matrix;

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  Vector<Real> * K_e = new Vector<Real>(0,
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
void Material::computePotentialEnergy(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if(ghost_type != _not_ghost) return;
  Real * epot = potential_energy(el_type, ghost_type)->values;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;

  computePotentialEnergy(strain_val, stress_val, epot);
  epot++;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::computePotentialEnergyByElement() {
  AKANTU_DEBUG_IN();

  const Mesh::ConnectivityTypeList & type_list = model->getFEM().getMesh().getConnectivityTypeList(_not_ghost);
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(model->getFEM().getMesh().getSpatialDimension(*it) != spatial_dimension) continue;

    if(!potential_energy.exists(*it, _not_ghost)) {
      UInt nb_element = element_filter(*it, _not_ghost)->getSize();
      UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(*it, _not_ghost);

      std::stringstream sstr; sstr << id << ":potential_energy:"<< *it;
      potential_energy(*it, _not_ghost) = &(alloc<Real> (sstr.str(), nb_element * nb_quadrature_points,
							 1,
							 REAL_INIT_VALUE));
    }

    computePotentialEnergy(*it);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real Material::getPotentialEnergy() {
  AKANTU_DEBUG_IN();
  Real epot = 0.;

  computePotentialEnergyByElement();

  /// integrate the potential energy for each type of elements
  const Mesh::ConnectivityTypeList & type_list = model->getFEM().getMesh().getConnectivityTypeList(_not_ghost);
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(model->getFEM().getMesh().getSpatialDimension(*it) != spatial_dimension) continue;

    epot += model->getFEM().integrate(*potential_energy(*it, _not_ghost), *it,
				      _not_ghost, element_filter(*it, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */
template void Material::initInternalVector<Real>(ByElementTypeVector<Real> & vect,
						 UInt nb_component,
						 const ID & id);
template void Material::initInternalVector<UInt>(ByElementTypeVector<UInt> & vect,
						 UInt nb_component,
						 const ID & id);
template void Material::initInternalVector<Int>(ByElementTypeVector<Int> & vect,
						UInt nb_component,
						const ID & id);

template void Material::resizeInternalVector<Real>(ByElementTypeVector<Real> & vect);
template void Material::resizeInternalVector<UInt>(ByElementTypeVector<UInt> & vect);
template void Material::resizeInternalVector<Int>(ByElementTypeVector<Int> & vect);

__END_AKANTU__
