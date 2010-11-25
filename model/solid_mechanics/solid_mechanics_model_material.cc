/**
 * @file   solid_mechanics_model_material.cc
 * @author Guillaume ANCIAUX <anciaux@lsmscluster1.epfl.ch>
 * @date   Thu Nov 25 10:48:53 2010
 *
 * @brief  
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "aka_math.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::readMaterials(const std::string & filename) {
  MaterialParser parser;
  parser.open(filename.c_str());
  std::string mat_type = parser.getNextMaterialType();
  UInt mat_count = 0;

  while (mat_type != ""){
    std::stringstream sstr_mat; sstr_mat << id << ":" << mat_count++ << ":" << mat_type;
    Material * material;
    MaterialID mat_id = sstr_mat.str();
    /// read the material properties
    if(mat_type == "elastic") material = parser.readMaterialObject<MaterialElastic>(*this,mat_id);
    else AKANTU_DEBUG_ERROR("Malformed material file : unknown material type "
			    << mat_type);
    materials.push_back(material);
    mat_type = parser.getNextMaterialType();
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initMaterials() {
  AKANTU_DEBUG_ASSERT(materials.size() != 0, "No material to initialize !");

  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    /// init internals properties
    (*mat_it)->initMaterial();
  }

  Material ** mat_val = &(materials.at(0));

  /// fill the element filters of the materials using the element_material arrays
  const Mesh::ConnectivityTypeList & type_list = fem->getMesh().getConnectivityTypeList(_not_ghost);
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    UInt nb_element = fem->getMesh().getNbElement(*it);
    UInt * elem_mat_val = element_material[*it]->values;

    for (UInt el = 0; el < nb_element; ++el) {
      mat_val[elem_mat_val[el]]->addElement(*it, el);
    }
  }

  /// @todo synchronize element material

  /// fill the element filters of the materials using the element_material arrays
  const Mesh::ConnectivityTypeList & ghost_type_list =
    fem->getMesh().getConnectivityTypeList(_ghost);
  for(it = ghost_type_list.begin(); it != ghost_type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    UInt nb_element = fem->getMesh().getNbGhostElement(*it);
    UInt * elem_mat_val = ghost_element_material[*it]->values;

    for (UInt el = 0; el < nb_element; ++el) {
      mat_val[elem_mat_val[el]]->addGhostElement(*it, el);
    }
  }
}

__END_AKANTU__
