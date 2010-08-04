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
inline void Material::setPotentialEnergyFlagOn(){
  AKANTU_DEBUG_IN();

  if(!potential_energy_vector) {
    /// for each connectivity types allocate the element filer array of the material
    UInt spatial_dimension = model->getSpatialDimension();
    const Mesh::ConnectivityTypeList & type_list = model->getFEM().getConnectivityTypeList();
    Mesh::ConnectivityTypeList::const_iterator it;
    for(it = type_list.begin(); it != type_list.end(); ++it) {
      if(model->getFEM().getSpatialDimension(*it) != spatial_dimension) continue;
      UInt nb_quadrature_points = model->getFEM().getNbQuadraturePoints(*it);
      UInt nb_element = element_filter[*it]->getSize();
      std::stringstream sstr; sstr << id << ":potential_energy:"<< *it;
      potential_energy[*it] = &(alloc<Real> (sstr.str(), nb_element, nb_quadrature_points, NAN));
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
inline const Vector<Real> & Material::getPotentialEnergy(ElementType type) const {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(potential_energy[type] != NULL,
		      "The material " << id << " has no element of kind : "<< type);
  AKANTU_DEBUG_OUT();
  return *potential_energy[type];
}
