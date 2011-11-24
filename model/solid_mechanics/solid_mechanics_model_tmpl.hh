/**
 * @file   solid_mechanics_model_tmpl.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Nov 24 08:45:57 2011
 *
 * @brief  template part of solid mechanics model
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
template <typename M>
UInt SolidMechanicsModel::readCustomMaterial(const std::string & filename,
					     const std::string & keyword) {

  Parser parser;
  parser.open(filename);
  std::string key = keyword;
  to_lower(key);
  std::string mat_name = parser.getNextSection("material");
  while (mat_name != ""){
    if (mat_name == key) break;
    mat_name = parser.getNextSection("material");
  }
  if (mat_name != key) AKANTU_DEBUG_ERROR("material "
					  << key
					  << " not found in file " << filename);

  std::stringstream sstr_mat; sstr_mat << id << ":" << materials.size() << ":" << key;
  ID mat_id = sstr_mat.str();
  Material * mat = parser.readSection<M>(*this, mat_id);
  materials.push_back(mat);
  return materials.size();;
}

/* -------------------------------------------------------------------------- */
