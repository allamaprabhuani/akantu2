/**
 * @file   constitutive_law.cc
 *
 * @author Mohit Pundir <mohit.pundir@ethz.ch>
 *
 * @date creation: Mon May 2 2022
 * @date last modification: Wed may 2 2022
 *
 * @brief  Implementation of common part of the constitutive law
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */


/* -------------------------------------------------------------------------- */
#include "constitutive_law.hh"
#include "poisson_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ConstitutiveLaw::ConstitutiveLaw(PoissonModel & model, const ID & id)
    : Parsable(ParserType::_constitutive_law, id), id(id), fem(model.getFEEngine()),
      model(model), spatial_dimension(this->model.getSpatialDimension()),
      element_filter("element_filter", id), 
      flux_dof("flux_degree_of_freedom", *this),
      gradient_dof("gradient_degree_of_freedom", *this)  {

  AKANTU_DEBUG_IN();

  /// for each connectivity types allocate the element filer array of the
  /// material
  element_filter.initialize(model.getMesh(),
                            _spatial_dimension = spatial_dimension,
                            _element_kind = _ek_regular);
  this->initialize();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
ConstitutiveLaw::ConstitutiveLaw(PoissonModel & model, UInt dim, const Mesh & mesh,
				 FEEngine & fe_engine, const ID & id)
    : Parsable(ParserType::_constitutive_law, id), id(id), fem(fe_engine),
      model(model), spatial_dimension(this->model.getSpatialDimension()),
      element_filter("element_filter", id),
      flux_dof("flux_degree_of_freedom", *this, dim, fe_engine,
	       this->element_filter),
      gradient_dof("gradient_degree_of_freedom", *this, dim, fe_engine,
		   this->element_filter) {

  AKANTU_DEBUG_IN();

  /// for each connectivity types allocate the element filer array of the
  /// material
  element_filter.initialize(mesh, _spatial_dimension = spatial_dimension,
                            _element_kind = _ek_regular);
  this->initialize();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ConstitutiveLaw::~ConstitutiveLaw() = default;


/* -------------------------------------------------------------------------- */
void ConstitutiveLaw::initialize() {
  registerParam("name", name, std::string(), _pat_parsable | _pat_readable);
 
  flux_dof.initialize(spatial_dimension);
  gradient_dof.initialize(spatial_dimension);
}

/* -------------------------------------------------------------------------- */
void ConstitutiveLaw::initConstitutiveLaw() {
  AKANTU_DEBUG_IN();

  this->resizeInternals();

  updateInternalParameters();

  AKANTU_DEBUG_OUT();
}


  
/* -------------------------------------------------------------------------- */
void ConstitutiveLaw::resizeInternals() {
  AKANTU_DEBUG_IN();
  for (auto it = internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it) {
    it->second->resize();
  }

  for (auto it = internal_vectors_uint.begin();
       it != internal_vectors_uint.end(); ++it) {
    it->second->resize();
  }

  for (auto it = internal_vectors_bool.begin();
       it != internal_vectors_bool.end(); ++it) {
    it->second->resize();
  }

  AKANTU_DEBUG_OUT();
}



/* -------------------------------------------------------------------------- */
void ConstitutiveLaw::computeAllFluxes(GhostType ghost_type) {

  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();

  for (const auto & type :
       element_filter.elementTypes(spatial_dimension, ghost_type)) {
    auto & elem_filter = element_filter(type, ghost_type);

    if (elem_filter.empty()) {
      continue;
    }

    computeFlux(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Constitutivelaw::assembleStiffnessMatrix(GhostType ghost_type) {

  AKANTU_DEBUG_IN();
  
  UInt spatial_dimension = model.getSpatialDimension();
  
  for (auto type : element_filter.elementTypes(spatial_dimension, ghost_type)) {

    Array<UInt> & elem_filter = element_filter(type, ghost_type);
    if (elem_filter.empty()) {
      AKANTU_DEBUG_OUT();
      return;
    }

    UInt nb_element = elem_filter.size();
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    UInt nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

    
    UInt tangent_size = spatial_dimension;
    
    auto * tangent_stiffness_matrix =
      new Array<Real>(nb_element * nb_quadrature_points,
                      tangent_size * tangent_size, "tangent_stiffness_matrix");
    
    tangent_stiffness_matrix->zero();
    
    computeTangentModuli(type, *tangent_stiffness_matrix, ghost_type);

    
    auto bt_d_b = std::make_unique<Array<Real>>(
				nb_element * nb_quadrature_points,
				nb_nodes_per_element * nb_nodes_per_element, "B^t*D*B");

    fem.computeBtDB(*tangent_stiffness_matrix, *bt_d_b, 2, type);

    delete tangent_stiffness_matrix;
    
    /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
    auto K_e = std::make_unique<Array<Real>>(
        nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_e");

    fem.integrate(*bt_d_b, *K_e, nb_nodes_per_element * nb_nodes_per_element,
                  type);

    delete bt_d_b;

    model.getDOFManager().assembleElementalMatricesToMatrix(
		  "K", "dof", *K_e, type, ghost_type, _symmetric, elem_filter);
    delete K_e;

  }

}
    
