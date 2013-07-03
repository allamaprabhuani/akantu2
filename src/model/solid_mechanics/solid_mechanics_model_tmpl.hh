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

/*
__END_AKANTU__
#include "sparse_matrix.hh"
#include "solver.hh"
__BEGIN_AKANTU__
*/
/* -------------------------------------------------------------------------- */
/*template<NewmarkBeta::IntegrationSchemeCorrectorType type>
void SolidMechanicsModel::solveDynamic(Array<Real> & increment) {
  AKANTU_DEBUG_INFO("Solving Ma + Cv + Ku = f");

  NewmarkBeta * nmb_int = dynamic_cast<NewmarkBeta *>(integrator);
  Real c = nmb_int->getAccelerationCoefficient<type>(time_step);
  Real d = nmb_int->getVelocityCoefficient<type>(time_step);
  Real e = nmb_int->getDisplacementCoefficient<type>(time_step);

  // A = c M + d C + e K
  jacobian_matrix->clear();

  if(stiffness_matrix)
    jacobian_matrix->add(*stiffness_matrix, e);

//  if(type != NewmarkBeta::_acceleration_corrector)
//    jacobian_matrix->add(*stiffness_matrix, e);

  if(mass_matrix)
    jacobian_matrix->add(*mass_matrix, c);

  mass_matrix->saveMatrix("M.mtx");
  if(velocity_damping_matrix)
    jacobian_matrix->add(*velocity_damping_matrix, d);

  jacobian_matrix->applyBoundary(*boundary);

#ifndef AKANTU_NDEBUG
  if(AKANTU_DEBUG_TEST(dblDump))
    jacobian_matrix->saveMatrix("J.mtx");
#endif

  jacobian_matrix->saveMatrix("J.mtx");

  solver->setRHS(*residual);

  // solve A w = f
  solver->solve(increment);
}*/
/* -------------------------------------------------------------------------- */
