/**
 * @file   material.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:43:41 2010
 *
 * @brief  Implementation of the common part of the material class
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
#include "material.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Material::Material(SolidMechanicsModel & model, const MaterialID & id) :
  Memory(model.getMemoryID()), id(id),
  spatial_dimension(model.getSpatialDimension()), name(""),
  model(&model), potential_energy_flag(false), potential_energy_vector(false),
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

  /// for each connectivity types allocate the element filer array of the material
  UInt spatial_dimension = model->getSpatialDimension();

  const Mesh::ConnectivityTypeList & type_list =
    model->getFEM().getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;
    std::stringstream sstr; sstr << id << ":element_filer:"<< *it;
    element_filter[*it] = &(alloc<UInt> (sstr.str(), 0, 1));

    std::stringstream sstr_stre; sstr_stre << id << ":stress:" << *it;
    std::stringstream sstr_stra; sstr_stra << id << ":strain:" << *it;

    stress[*it] = &(alloc<Real>(sstr_stre.str(), 0,
				spatial_dimension * spatial_dimension));
    strain[*it] = &(alloc<Real>(sstr_stra.str(), 0,
				spatial_dimension * spatial_dimension));
  }


  const Mesh::ConnectivityTypeList & ghost_type_list =
    model->getFEM().getMesh().getConnectivityTypeList(_ghost);

  for(it = ghost_type_list.begin(); it != ghost_type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    std::stringstream sstr; sstr << id << ":ghost_element_filer:"<< *it;
    ghost_element_filter[*it] = &(alloc<UInt> (sstr.str(), 0, 1));

    std::stringstream sstr_stre; sstr_stre << id << ":ghost_stress:" << *it;
    std::stringstream sstr_stra; sstr_stra << id << ":ghost_strain:" << *it;

    ghost_stress[*it] = &(alloc<Real>(sstr_stre.str(), 0,
				spatial_dimension * spatial_dimension));
    ghost_strain[*it] = &(alloc<Real>(sstr_stra.str(), 0,
				spatial_dimension * spatial_dimension));
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
    UInt size_of_shapes_derivatives = FEM::getShapeDerivativesSize(*it);
    UInt nb_quadrature_points       = FEM::getNbQuadraturePoints(*it);

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

    UInt nb_element = elem_filter->getSize();

    /// compute @f$\nabla u@f$
    model->getFEM().gradientOnQuadraturePoints(current_position, *strain_vect,
					      spatial_dimension,
					      *it, ghost_type, elem_filter);

    /// compute @f$\mathbf{\sigma}_q@f$ from @f$\nabla u@f$
    computeStress(*it, ghost_type);

    if(potential_energy_flag) computePotentialEnergy(*it, ghost_type);

    /// compute @f$\sigma \frac{\partial \varphi}{\partial X}@f$ by @f$\mathbf{B}^t \mathbf{\sigma}_q@f$
    Vector<Real> * sigma_dphi_dx =
      new Vector<Real>(nb_element*nb_quadrature_points, size_of_shapes_derivatives, "sigma_x_dphi_/_dX");

    Real * shapesd           = shapes_derivatives->values;
    UInt size_of_shapesd     = shapes_derivatives->getNbComponent();
    Real * shapesd_val;
    Real * stress_val        = stress_vect->values;
    Real * sigma_dphi_dx_val = sigma_dphi_dx->values;
    UInt * elem_filter_val = elem_filter->values;

    UInt offset_shapesd_val       = spatial_dimension * nb_nodes_per_element;
    UInt offset_stress_val        = spatial_dimension * spatial_dimension;
    UInt offset_sigma_dphi_dx_val = spatial_dimension * nb_nodes_per_element;

    for (UInt el = 0; el < nb_element; ++el) {
      shapesd_val = shapesd + elem_filter_val[el]*size_of_shapesd*nb_quadrature_points;
      for (UInt q = 0; q < nb_quadrature_points; ++q) {
	Math::matrix_matrixt(nb_nodes_per_element, spatial_dimension, spatial_dimension,
			     shapesd_val, stress_val, sigma_dphi_dx_val);
	shapesd_val       += offset_shapesd_val;
	stress_val        += offset_stress_val;
	sigma_dphi_dx_val += offset_sigma_dphi_dx_val;
      }
    }

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
Real Material::getPotentialEnergy() {
  AKANTU_DEBUG_IN();
  Real epot = 0.;

  /// integrate the potential energy for each type of elements
  if(potential_energy_flag) {
    const Mesh::ConnectivityTypeList & type_list = model->getFEM().getMesh().getConnectivityTypeList(_not_ghost);
    Mesh::ConnectivityTypeList::const_iterator it;
    for(it = type_list.begin(); it != type_list.end(); ++it) {
      if(model->getFEM().getMesh().getSpatialDimension(*it) != spatial_dimension) continue;

      epot += model->getFEM().integrate(*potential_energy[*it], *it,
					_not_ghost, element_filter[*it]);
    }
  }

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
