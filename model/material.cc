/**
 * @file   material.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:43:41 2010
 *
 * @brief  Implementation of the common part of the material class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "material.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Material::Material(SolidMechanicsModel & model, const MaterialID & id) :
  Memory(model.getMemoryID()), id(id), model(&model),
  potential_energy_flag(false), potential_energy_vector(false),
  is_init(false) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Material::~Material() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model->getSpatialDimension();

  const Mesh::ConnectivityTypeList & type_list = model->getFEM().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(model->getFEM().getSpatialDimension(*it) != spatial_dimension) continue;
    dealloc(element_filter[*it]->getID());
    element_filter[*it] = NULL;

    if(potential_energy[*it]){
      dealloc(potential_energy[*it]->getID());
      potential_energy[*it] = NULL;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::initMaterial() {
  AKANTU_DEBUG_IN();

  /// for each connectivity types allocate the element filer array of the material
  UInt spatial_dimension = model->getSpatialDimension();

  const Mesh::ConnectivityTypeList & type_list = model->getFEM().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(model->getFEM().getSpatialDimension(*it) != spatial_dimension) continue;
    std::stringstream sstr; sstr << id << ":element_filer:"<< *it;
    element_filter[*it] = &(alloc<UInt> (sstr.str(), 0, 1));
  }

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
