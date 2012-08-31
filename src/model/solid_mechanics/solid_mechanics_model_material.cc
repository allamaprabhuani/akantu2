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

/* -------------------------------------------------------------------------- */
#define AKANTU_INTANTIATE_MATERIAL_BY_DIM_NO_TMPL(dim, elem)		\
  material =								\
    &(registerNewCustomMaterial< BOOST_PP_ARRAY_ELEM(1, elem)< dim > >(mat_type, \
								      opt_param))

#define AKANTU_INTANTIATE_MATERIAL_BY_DIM_TMPL_EACH(r, data, i, elem)	\
  BOOST_PP_EXPR_IF(BOOST_PP_NOT_EQUAL(0, i), else )			\
    if(opt_param == BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(2, 0, elem))) { \
      material =							\
    &(registerNewCustomMaterial< BOOST_PP_ARRAY_ELEM(1, data)< BOOST_PP_ARRAY_ELEM(0, data), \
							       BOOST_PP_TUPLE_ELEM(2, 1, elem) > >(mat_type, opt_param)); \
    }


#define AKANTU_INTANTIATE_MATERIAL_BY_DIM_TMPL(dim, elem)		\
  BOOST_PP_SEQ_FOR_EACH_I(AKANTU_INTANTIATE_MATERIAL_BY_DIM_TMPL_EACH,	\
			  (2, (dim, BOOST_PP_ARRAY_ELEM(1, elem))),	\
			  BOOST_PP_ARRAY_ELEM(2, elem))			\
  else {								\
    AKANTU_INTANTIATE_MATERIAL_BY_DIM_NO_TMPL(dim, elem);		\
  }

#define AKANTU_INTANTIATE_MATERIAL_BY_DIM(dim, elem)			\
  BOOST_PP_IF(BOOST_PP_EQUAL(3, BOOST_PP_ARRAY_SIZE(elem) ),		\
   	      AKANTU_INTANTIATE_MATERIAL_BY_DIM_TMPL,			\
    	      AKANTU_INTANTIATE_MATERIAL_BY_DIM_NO_TMPL)(dim, elem)


#define AKANTU_INTANTIATE_MATERIAL(elem)				\
  switch(spatial_dimension) {						\
  case 1: { AKANTU_INTANTIATE_MATERIAL_BY_DIM(1, elem); break; }	\
  case 2: { AKANTU_INTANTIATE_MATERIAL_BY_DIM(2, elem); break; }	\
  case 3: { AKANTU_INTANTIATE_MATERIAL_BY_DIM(3, elem); break; }	\
  }


#define AKANTU_INTANTIATE_MATERIAL_IF(elem)			\
  if (mat_type == BOOST_PP_STRINGIZE(BOOST_PP_ARRAY_ELEM(0, elem))) { \
    AKANTU_INTANTIATE_MATERIAL(elem);					\
  }

#define AKANTU_INTANTIATE_OTHER_MATERIAL(r, data, elem)			\
  else AKANTU_INTANTIATE_MATERIAL_IF(elem)

#define AKANTU_INSTANTIATE_MATERIALS()					\
  do {									\
    AKANTU_INTANTIATE_MATERIAL_IF(BOOST_PP_SEQ_HEAD(AKANTU_MATERIAL_LIST)) \
      BOOST_PP_SEQ_FOR_EACH(AKANTU_INTANTIATE_OTHER_MATERIAL,		\
			    _,						\
			    BOOST_PP_SEQ_TAIL(AKANTU_MATERIAL_LIST))	\
    else								\
      AKANTU_DEBUG_ERROR("Malformed material file : unknown material type " \
			 << mat_type);					\
  } while(0)

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::readMaterials(const std::string & filename) {
  Parser parser;
  parser.open(filename);
  std::string opt_param;
  std::string mat_type = parser.getNextSection("material", opt_param);

  while (mat_type != ""){
    Material & mat = registerNewMaterial(mat_type, opt_param);
    parser.readSection(mat.getID(), mat);
    mat_type = parser.getNextSection("material", opt_param);
  }
  parser.close();
}

/* -------------------------------------------------------------------------- */
Material & SolidMechanicsModel::registerNewMaterial(const ID & mat_type,
						    const std::string & opt_param) {
  // UInt mat_count = materials.size();

  // std::stringstream sstr_mat; sstr_mat << id << ":" << mat_count << ":" << mat_type;
  // Material * material;
  // ID mat_id = sstr_mat.str();

  // // add all the new materials in the AKANTU_MATERIAL_LIST in the material.hh file
  // AKANTU_INSTANTIATE_MATERIALS();
  // materials.push_back(material);
  Material * material = NULL;

  AKANTU_INSTANTIATE_MATERIALS();

  return *material;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initMaterials() {
  AKANTU_DEBUG_ASSERT(materials.size() != 0, "No material to initialize !");

  Material ** mat_val = &(materials.at(0));

  /// @todo synchronize element material

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    /// fill the element filters of the materials using the element_material arrays
    Mesh::type_iterator it  = mesh.firstType(spatial_dimension, gt);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt);

    for(; it != end; ++it) {
      UInt nb_element = mesh.getNbElement(*it, gt);
      UInt * elem_mat_val = element_material(*it, gt).storage();
      element_index_by_material.alloc(nb_element, 1, *it, gt);

      for (UInt el = 0; el < nb_element; ++el) {
	UInt index = mat_val[elem_mat_val[el]]->addElement(*it, el, gt);
	element_index_by_material(*it, gt)(el) = index;
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


/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::setMaterialIDsFromIntData(const std::string & data_name) {
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    Mesh::type_iterator it  = mesh.firstType(spatial_dimension, gt);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt);

    for(; it != end; ++it) {
      try {
	const Vector<UInt> & data = mesh.getUIntData(*it, data_name, gt);

	AKANTU_DEBUG_ASSERT(element_material.exists(*it, gt),
			    "element_material for type (" << gt << ":" << *it
			    << ") does not exists!");

	element_material(*it, gt).copy(data);
      } catch(...) {
	AKANTU_DEBUG_ERROR("No data named " << data_name
			   << " present in the mesh " << id
			   << " for the element type " << *it);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
