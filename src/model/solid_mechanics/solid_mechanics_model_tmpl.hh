/**
 * @file   solid_mechanics_model_tmpl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 *
 * @date   Thu Nov 24 09:36:33 2011
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
Material & SolidMechanicsModel::registerNewCustomMaterial(const ID & mat_type,
							  __attribute__((unused)) const std::string & opt_param) {
  UInt mat_count = materials.size();

  std::stringstream sstr_mat; sstr_mat << id << ":" << mat_count << ":" << mat_type;
  Material * material;
  ID mat_id = sstr_mat.str();

  // add all the new materials in the AKANTU_MATERIAL_LIST in the material.hh file
  material = new M(*this, mat_id);
  materials.push_back(material);

  return *material;
}



/* -------------------------------------------------------------------------- */
template <typename M>
UInt SolidMechanicsModel::readCustomMaterial(const std::string & filename,
					     const std::string & keyword) {

  Parser parser;
  parser.open(filename);
  std::string key = keyword;

  std::string opt_param;
  std::string mat_name = parser.getNextSection("material", opt_param);
  while (mat_name != ""){
    if (mat_name == key) break;
    mat_name = parser.getNextSection("material", opt_param);
  }

  if (mat_name != key) AKANTU_DEBUG_ERROR("material "
					  << key
					  << " not found in file " << filename);

  Material & mat = registerNewCustomMaterial<M>(key, opt_param);
  parser.readSection(mat.getID(), mat);
  materials.push_back(&mat);
  return materials.size();;
}

/* --------------------------------------------------------------------------*/
//template<class Functor>
//void SolidMechanicsModel::computeForcesFromFunction(Functor & functor,
						    //BoundaryFunctionType function_type) {
  ///** function type is
   //** _bft_forces : traction function is given
   //** _bft_stress : stress function is given
   //*/
  //GhostType ghost_type = _not_ghost;

  //UInt nb_component = 0;
  //switch(function_type) {
  //case _bft_stress: nb_component = spatial_dimension * spatial_dimension; break;
  //case _bft_traction: nb_component = spatial_dimension; break;
  //default: break;
  //}

  //Array<Real> funct(0, nb_component, "traction_stress");
  //Array<Real> quad_coords(0, spatial_dimension, "quad_coords");

  ////prepare the loop over element types
  //Mesh::type_iterator it  = getFEMBoundary().getMesh().firstType(getFEMBoundary().getElementDimension(),
								 //ghost_type);
  //Mesh::type_iterator end = getFEMBoundary().getMesh().lastType(getFEMBoundary().getElementDimension(),
								//ghost_type);

  //for(; it != end; ++it) {

    //UInt nb_quad    = getFEMBoundary().getNbQuadraturePoints(*it, ghost_type);
    //UInt nb_element = getFEMBoundary().getMesh().getNbElement(*it, ghost_type);

    //funct.resize(nb_element * nb_quad);
    //quad_coords.resize(nb_element * nb_quad);

    //const Array<Real> & normals_on_quad = getFEMBoundary().getNormalsOnQuadPoints(*it, ghost_type);

    //getFEMBoundary().interpolateOnQuadraturePoints(getFEMBoundary().getMesh().getNodes(),
						   //quad_coords, spatial_dimension, *it, ghost_type);

    //Array<Real>::const_iterator< Vector<Real> > normals = normals_on_quad.begin(spatial_dimension);
    //Array<Real>::iterator< Vector<Real> > qcoord  = quad_coords.begin(spatial_dimension);


    //Array<UInt>::iterator< UInt > surface_id;
    //bool has_surface_id;
    //// XXX TODO FIXME
////    try {
////      surface_id = mesh.getSurfaceID(*it, ghost_type).begin();
////      has_surface_id = true;
////    } catch (...) {
////      has_surface_id = false;
////    }
////
    //if(function_type == _bft_stress) {
      //Array<Real>::iterator< Matrix<Real> > stress = funct.begin(spatial_dimension, spatial_dimension);

      //for (UInt el = 0; el < nb_element; ++el) {
	      //Surface surf_id = 0;
	      //if(has_surface_id) {
	        //surf_id = *surface_id;
	        //++surface_id;
	      //}
	      //for (UInt q = 0; q < nb_quad; ++q, ++stress, ++qcoord, ++normals) {
	        //functor.stress(*qcoord, *stress, *normals, surf_id);
	      //}
      //}
    //} else if (function_type == _bft_traction) {
      //Array<Real>::iterator< Vector<Real> > force = funct.begin(spatial_dimension);

      //for (UInt el = 0; el < nb_element; ++el) {
	//Surface surf_id = 0;
	//if(has_surface_id) {
	  //surf_id = *surface_id;
	  //++surface_id;
	//}
	//for (UInt q = 0; q < nb_quad; ++q, ++force, ++qcoord, ++normals) {
	  //functor.traction(*qcoord, *force, *normals, surf_id);
	//}
      //}
    //}

    //switch(function_type) {
    //case _bft_stress:
      //computeForcesByStressTensor(funct, *it, ghost_type); break;
    //case _bft_traction:
      //computeForcesByTractionArray(funct, *it, ghost_type); break;
    //default: break;
    //}
  //}
//}

/* -------------------------------------------------------------------------- */

