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

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Material::Material(SolidMechanicsModel & model, const MaterialID & id) :
  Memory(model.getMemoryID()), id(id),
  spatial_dimension(model.getSpatialDimension()), name(""),
  model(&model), potential_energy_vector(false),
  is_init(false) {
  AKANTU_DEBUG_IN();

  for(UInt t = _not_defined; t < _max_element_type; ++t) {
    this->stress                [t] = NULL;
    this->strain                [t] = NULL;
    this->potential_energy      [t] = NULL;
    this->element_filter        [t] = NULL;
    this->ghost_stress          [t] = NULL;
    this->ghost_strain          [t] = NULL;
    this->ghost_potential_energy[t] = NULL;
    this->ghost_element_filter  [t] = NULL;
  }

  /// allocate strain stress for local elements
  initInternalVector(strain,       spatial_dimension*spatial_dimension, "strain", _not_ghost);
  initInternalVector(stress,       spatial_dimension*spatial_dimension, "stress", _not_ghost);

  /// allocate strain stress for ghost elements
  initInternalVector(ghost_strain, spatial_dimension*spatial_dimension, "strain", _ghost);
  initInternalVector(ghost_stress, spatial_dimension*spatial_dimension, "stress", _ghost);

  /// for each connectivity types allocate the element filer array of the material
  const Mesh::ConnectivityTypeList & type_list =
    model.getFEM().getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;
    std::stringstream sstr; sstr << id << ":element_filer:"<< *it;
    element_filter[*it] = &(alloc<UInt> (sstr.str(), 0, 1));
  }


  const Mesh::ConnectivityTypeList & ghost_type_list =
    model.getFEM().getMesh().getConnectivityTypeList(_ghost);

  for(it = ghost_type_list.begin(); it != ghost_type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    std::stringstream sstr; sstr << id << ":ghost_element_filer:"<< *it;
    ghost_element_filter[*it] = &(alloc<UInt> (sstr.str(), 0, 1));
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
  resizeInternalVector(stress,_not_ghost);
  resizeInternalVector(ghost_stress,_ghost);
  resizeInternalVector(strain,_not_ghost);
  resizeInternalVector(ghost_strain,_ghost);
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void Material::initInternalVector(ByElementTypeReal & vect,
				  UInt nb_component,
				  const std::string & vect_id,
				  GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  model->getFEM().getMesh().initByElementTypeRealVector(vect,nb_component,spatial_dimension,
							 id,vect_id,ghost_type);
  // for(UInt t = _not_defined; t < _max_element_type; ++t)
  //   vect[t] = NULL;

  // std::string ghost_id = "";

  // if (ghost_type == _ghost) {
  //   ghost_id = "ghost_";
  // }

  // const Mesh::ConnectivityTypeList & type_list = model->getFEM().getMesh().getConnectivityTypeList();
  // Mesh::ConnectivityTypeList::const_iterator it;
  // for(it = type_list.begin(); it != type_list.end(); ++it) {
  //   if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;
  //   std::stringstream sstr_damage; sstr_damage << id << ":" << ghost_id << vect_id << ":" << *it;
  //   vect[*it] = &(alloc<Real>(sstr_damage.str(), 0,
  // 			      nb_component, REAL_INIT_VALUE));
  // }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void Material::resizeInternalVector(ByElementTypeReal & vect,
				    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const Mesh::ConnectivityTypeList & type_list = model->getFEM().getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    Vector<UInt> * elem_filter;
    if (ghost_type == _not_ghost) {
      elem_filter = element_filter[*it];
    } else if (ghost_type == _ghost) {
      elem_filter = ghost_element_filter[*it];
    }

    UInt nb_element           = elem_filter->getSize();
    UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(*it);
    UInt new_size = nb_element*nb_quadrature_points;

    UInt size = vect[*it]->getSize();
    UInt nb_component = vect[*it]->getNbComponent();
    vect[*it]->resize(new_size);
    memset(vect[*it]->values + size * nb_component, 0, (new_size - size) * nb_component * sizeof(Real));
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

    UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(*it);
    UInt nb_quadrature_points       = model->getFEM().getNbQuadraturePoints(*it);

    Vector<Real> * strain_vect;
    Vector<Real> * stress_vect;
    Vector<UInt> * elem_filter;
    const Vector<Real> * shapes_derivatives;

    if(ghost_type == _not_ghost) {
      elem_filter = element_filter[*it];
      strain_vect = strain[*it];
      stress_vect = stress[*it];
      shapes_derivatives = &(model->getFEM().getShapesDerivatives(*it));
    } else {
      elem_filter = ghost_element_filter[*it];
      strain_vect = ghost_strain[*it];
      stress_vect = ghost_stress[*it];
      shapes_derivatives = &(model->getFEM().getGhostShapesDerivatives(*it));
    }

    UInt size_of_shapes_derivatives = shapes_derivatives->getNbComponent();

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

    Real * shapesd           = shapes_derivatives->values;
    UInt size_of_shapesd     = shapes_derivatives->getNbComponent();
    Real * shapesd_val;
    //    Real * stress_val        = stress_vect->values;
    //    Real * sigma_dphi_dx_val = sigma_dphi_dx->values;
    UInt * elem_filter_val   = elem_filter->values;

    //    UInt offset_shapesd_val       = spatial_dimension * nb_nodes_per_element;
    //    UInt offset_stress_val        = spatial_dimension * spatial_dimension;
    //    UInt offset_sigma_dphi_dx_val = spatial_dimension * nb_nodes_per_element;

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
    Vector<Real> * int_sigma_dphi_dx = new Vector<Real>(0,
							nb_nodes_per_element * spatial_dimension,
							"int_sigma_x_dphi_/_dX");
    model->getFEM().integrate(*sigma_dphi_dx, *int_sigma_dphi_dx,
			      size_of_shapes_derivatives,
			      *it, ghost_type,
			      elem_filter);
    delete sigma_dphi_dx;

    /// assemble
    model->getFEM().assembleVector(*int_sigma_dphi_dx, residual,
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
void Material::computeStiffnessMatrix(Vector<Real> & current_position, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model->getSpatialDimension();

  SparseMatrix K = const_cast<SparseMatrix &>(model->getStiffnessMatrix());

  const Mesh::ConnectivityTypeList & type_list =
    model->getFEM().getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {

    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(*it);
    UInt nb_quadrature_points       = model->getFEM().getNbQuadraturePoints(*it);

    Vector<Real> * strain_vect;
    Vector<UInt> * elem_filter;
    const Vector<Real> * shapes_derivatives;

    if(ghost_type == _not_ghost) {
      elem_filter = element_filter[*it];
      strain_vect = strain[*it];
      shapes_derivatives = &(model->getFEM().getShapesDerivatives(*it));
    } else {
      elem_filter = ghost_element_filter[*it];
      strain_vect = ghost_strain[*it];
      shapes_derivatives = &(model->getFEM().getGhostShapesDerivatives(*it));
    }

    UInt size_of_shapes_derivatives = shapes_derivatives->getNbComponent();

    UInt nb_element = elem_filter->getSize();

    model->getFEM().gradientOnQuadraturePoints(current_position, *strain_vect,
					      spatial_dimension,
					      *it, ghost_type, elem_filter);

    //computeTangentStiffness(*it, ghost_type);

    /// compute @f$\sigma \frac{\partial \varphi}{\partial X}@f$ by @f$\mathbf{B}^t \mathbf{\sigma}_q@f$
    Vector<Real> * bt_d_b =
      new Vector<Real>(nb_element*nb_quadrature_points, size_of_shapes_derivatives, "sigma_x_dphi_/_dX");

    Real * shapesd           = shapes_derivatives->values;
    UInt size_of_shapesd     = shapes_derivatives->getNbComponent();
    // Real * shapesd_val;
    //Real * stress_val        = stress_vect->values;
    // Real * sigma_dphi_dx_val = sigma_dphi_dx->values;
    // UInt * elem_filter_val   = elem_filter->values;

    

    // UInt offset_shapesd_val       = spatial_dimension * nb_nodes_per_element;
    // UInt offset_stress_val        = spatial_dimension * spatial_dimension;
    // UInt offset_sigma_dphi_dx_val = spatial_dimension * nb_nodes_per_element;

    // Vector<Real> * shapesd_filtered =
    //   new Vector<Real>(nb_element*nb_quadrature_points, size_of_shapes_derivatives, "filtered shapesd");
    // Real * shapesd_filtered_val = shapesd_filtered->values;

    // for (UInt el = 0; el < nb_element; ++el) {
    //   shapesd_val = shapesd + elem_filter_val[el] * size_of_shapesd * nb_quadrature_points;
    //   memcpy(shapesd_filtered_val,
    // 	     shapesd_val,
    // 	     offset_shapesd_val * nb_quadrature_points * sizeof(Real));
    //   shapesd_filtered_val += offset_shapesd_val * nb_quadrature_points;
    // }

    // Math::matrix_matrixt(nb_nodes_per_element, spatial_dimension, spatial_dimension,
    // 			 *shapesd_filtered,
    // 			 *stress_vect,
    // 			 *sigma_dphi_dx);

    // delete shapesd_filtered;

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
    // Vector<Real> * int_sigma_dphi_dx = new Vector<Real>(0,
    // 							nb_nodes_per_element * spatial_dimension,
    // 							"int_sigma_x_dphi_/_dX");
    // model->getFEM().integrate(*sigma_dphi_dx, *int_sigma_dphi_dx,
    // 			      size_of_shapes_derivatives,
    // 			      *it, ghost_type,
    // 			      elem_filter);
    // delete sigma_dphi_dx;

    /// assemble
    // model->getFEM().assembleVector(*int_sigma_dphi_dx, residual,
    // 				   residual.getNbComponent(),
    // 				   *it, ghost_type, elem_filter, -1);
    //    delete int_sigma_dphi_dx;
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void Material::computePotentialEnergyByElement() {
  AKANTU_DEBUG_IN();

  const Mesh::ConnectivityTypeList & type_list = model->getFEM().getMesh().getConnectivityTypeList(_not_ghost);
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(model->getFEM().getMesh().getSpatialDimension(*it) != spatial_dimension) continue;

    if(potential_energy[*it] == NULL) {
      UInt nb_element = element_filter[*it]->getSize();
      UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(*it);

      std::stringstream sstr; sstr << id << ":potential_energy:"<< *it;
      potential_energy[*it] = &(alloc<Real> (sstr.str(), nb_element * nb_quadrature_points,
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

    epot += model->getFEM().integrate(*potential_energy[*it], *it,
				      _not_ghost, element_filter[*it]);
  }

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
