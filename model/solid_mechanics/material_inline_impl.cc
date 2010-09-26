/**
 * @file   material_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:57:43 2010
 *
 * @brief  Implementation of the inline functions of the class material
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */


/* -------------------------------------------------------------------------- */
inline void Material::addElement(ElementType type, UInt element) {
  element_filter[type]->push_back(element);
  strain[type]->push_back(REAL_INIT_VALUE);
  stress[type]->push_back(REAL_INIT_VALUE);
  if(potential_energy_vector)
    potential_energy[type]->push_back(REAL_INIT_VALUE);
}

/* -------------------------------------------------------------------------- */
inline void Material::addGhostElement(ElementType type, UInt element) {
  ghost_element_filter[type]->push_back(element);
  ghost_strain[type]->push_back(REAL_INIT_VALUE);
  ghost_stress[type]->push_back(REAL_INIT_VALUE);

  // if(potential_energy_vector)
  //   ghost_potential_energy[type]->push_back(REAL_INIT_VALUE);
}

/* -------------------------------------------------------------------------- */
inline void Material::setPotentialEnergyFlagOn(){
  AKANTU_DEBUG_IN();

  if(!potential_energy_vector) {
    /// for each connectivity types allocate the element filer array of the material
    for(UInt t = _not_defined + 1; t < _max_element_type; ++t) {
      ElementType type = (ElementType) t;
      if(Mesh::getSpatialDimension(type) != spatial_dimension) continue;

      UInt nb_quadrature_points = FEM::getNbQuadraturePoints(type);
      if(element_filter[type] != NULL) {
	UInt nb_element = element_filter[type]->getSize();
	std::stringstream sstr; sstr << id << ":potential_energy:"<< type;
	potential_energy[type] = &(alloc<Real> (sstr.str(), nb_element,
						nb_quadrature_points,
						REAL_INIT_VALUE));
      }

      // if(ghost_element_filter[type] != NULL) {
      // 	UInt nb_element = ghost_element_filter[type]->getSize();
      // 	std::stringstream sstr; sstr << id << ":ghost_potential_energy:"<< type;
      // 	ghost_potential_energy[type] = &(alloc<Real> (sstr.str(), nb_element,
      // 						      nb_quadrature_points,
      // 						      REAL_INIT_VALUE));
      // }
    }
    potential_energy_vector = true;
  }

  potential_energy_flag = true;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void Material::setPotentialEnergyFlagOff(){
  AKANTU_DEBUG_IN();
  potential_energy_flag = false;
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
inline UInt Material::getNbDataToPack(const Element & element,
				      GhostSynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
  return 0;
};

/* -------------------------------------------------------------------------- */
inline UInt Material::getNbDataToUnpack(const Element & element,
					GhostSynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
  return 0;
};

/* -------------------------------------------------------------------------- */
inline void Material::packData(Real ** buffer,
			       const Element & element,
			       GhostSynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
};

/* -------------------------------------------------------------------------- */
inline void Material::unpackData(Real ** buffer,
				 const Element & element,
				 GhostSynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
};
