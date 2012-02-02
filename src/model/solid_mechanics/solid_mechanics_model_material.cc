/**
 * @file   solid_mechanics_model_material.cc
 * @author Guillaume ANCIAUX <guillaume.anciaux@epfl.ch>
 * @date   Thu Nov 25 10:48:53 2010
 *
 * @brief
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
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "aka_math.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
#define AKANTU_INTANTIATE_OTHER_MATERIAL(r, data, elem)			\
  else if (BOOST_PP_TUPLE_ELEM(4, 0, data) ==				\
	   BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(2, 0, elem)))		\
    BOOST_PP_TUPLE_ELEM(4, 1, data) =					\
      BOOST_PP_TUPLE_ELEM(4, 3, data).					\
      readSection<BOOST_PP_TUPLE_ELEM(2, 1, elem)>(*this,		\
						   BOOST_PP_TUPLE_ELEM(4, 2, data));

#define AKANTU_INTANTIATE_OTHER_MATERIALS(mat_type, material, mat_id, parser, mat_lst) \
  BOOST_PP_SEQ_FOR_EACH(AKANTU_INTANTIATE_OTHER_MATERIAL,		\
			(mat_type, material, mat_id, parser),		\
			mat_lst)

#define AKANTU_INTANTIATE_MATERIALS(mat_type, material, mat_id, parser)	\
  do {									\
    if(mat_type ==							\
       BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(2,			\
					      0,			\
					      BOOST_PP_SEQ_HEAD(AKANTU_MATERIAL_LIST)))) \
      material =							\
	parser.readSection<BOOST_PP_TUPLE_ELEM(2,			\
					       1,			\
					       BOOST_PP_SEQ_HEAD(AKANTU_MATERIAL_LIST))> \
	(*this, mat_id);						\
    AKANTU_INTANTIATE_OTHER_MATERIALS(mat_type, material, mat_id, parser, \
				      BOOST_PP_SEQ_TAIL(AKANTU_MATERIAL_LIST)) \
    else AKANTU_DEBUG_ERROR("Malformed material file : unknown material type " \
			    << mat_type);				\
  } while(0)
/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::readMaterials(const std::string & filename) {
  Parser parser;
  parser.open(filename);
  std::string mat_type = parser.getNextSection("material");
  UInt mat_count = 0;

  while (mat_type != ""){
    std::stringstream sstr_mat; sstr_mat << id << ":" << mat_count++ << ":" << mat_type;
    Material * material;
    ID mat_id = sstr_mat.str();
    /// read the material properties

    // add all the new materials in the AKANTU_MATERIAL_LIST in the material.hh file
    AKANTU_INTANTIATE_MATERIALS(mat_type, material, mat_id, parser);

    materials.push_back(material);
    mat_type = parser.getNextSection("material");
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initMaterials() {
  AKANTU_DEBUG_ASSERT(materials.size() != 0, "No material to initialize !");

  Material ** mat_val = &(materials.at(0));

  /// @todo synchronize element material

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    /// fill the element filters of the materials using the element_material arrays
    const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList(gt);
    Mesh::ConnectivityTypeList::const_iterator it;
    for(it = type_list.begin(); it != type_list.end(); ++it) {
      if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

      UInt nb_element = mesh.getNbElement(*it, gt);
      UInt * elem_mat_val = element_material(*it, gt).storage();

      for (UInt el = 0; el < nb_element; ++el) {
	mat_val[elem_mat_val[el]]->addElement(*it, el, gt);
      }
    }
  }

  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    /// init internals properties
    (*mat_it)->initMaterial();
  }

  synch_registry->synchronize(_gst_smm_init_mat);
}

__END_AKANTU__
