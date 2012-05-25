/**
 * @file   material_neohookean.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:53:52 2010
 *
 * @brief  Specialization of the material class for the elastic material
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
#include "material_neohookean.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialNeohookean<spatial_dimension>::MaterialNeohookean(SolidMechanicsModel & model,
							  const ID & id)  :
  Material(model, id) {
  AKANTU_DEBUG_IN();

  rho = 0;
  E   = 0;
  nu  = 1./2.;
  plane_stress = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNeohookean<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();

  lambda   = nu * E / ((1 + nu) * (1 - 2*nu));
  mu       = E / (2 * (1 + nu));
  if(plane_stress)
    lambda = 2 * lambda * mu / (lambda + 2 * mu);
  kpa      = lambda + 2./3. * mu;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNeohookean<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real F[3*3];
  Real sigma[3*3];

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;
  memset(F, 0, 3 * 3 * sizeof(Real));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      F[3*i + j] = strain_val[spatial_dimension * i + j];
  for (UInt i = 0; i < spatial_dimension; ++i) F[i*3 + i] += 1;
  //F is now the real F !!

  computeStress(F, sigma);

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      stress_val[spatial_dimension*i + j] = sigma[3 * i + j];

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNeohookean<spatial_dimension>::computePotentialEnergy(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if(ghost_type != _not_ghost) return;
  Real * epot = potential_energy(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;

  computePotentialEnergy(strain_val, epot);
  epot++;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNeohookean<spatial_dimension>::computeTangentStiffness(const ElementType & el_type,
								    Vector<Real> & tangent_matrix,
								    GhostType ghost_type) {
//   switch(spatial_dimension) {
//   case 1: { computeTangentStiffnessByDim<1>(el_type, tangent_matrix, ghost_type); break; }
//   case 2: { computeTangentStiffnessByDim<2>(el_type, tangent_matrix, ghost_type); break; }
//   case 3: { computeTangentStiffnessByDim<3>(el_type, tangent_matrix, ghost_type); break; }
//   }
// }


// /* -------------------------------------------------------------------------- */
// template<UInt dim>
// void  MaterialNeohookean::computeTangentStiffnessByDim(ElementType el_type,
// 						       Vector<Real>& tangent_matrix,
// 						       GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  tangent_matrix.clear();
  Real F[3*3];
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);

  memset(F, 0, 3 * 3 * sizeof(Real));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      F[3*i + j] = strain_val[spatial_dimension * i + j];
  for (UInt i = 0; i < spatial_dimension; ++i) F[i*3 + i] += 1;
  //F is now the real F !!

  computeTangentStiffness(tangent_val, F);
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
Real MaterialNeohookean<spatial_dimension>::celerity(const Element & elem) {
  UInt nb_quadrature_points =
    model->getFEM().getNbQuadraturePoints(elem.type, elem.ghost_type);

  UInt size_strain = spatial_dimension * spatial_dimension;

  Real * strain_val = strain(elem.type, elem.ghost_type).storage()
    + elem.element * nb_quadrature_points * size_strain;

  Real F[3*3];
  Real C[3*3], Cinv[3*3];

  Real cele = 0.;

  for (UInt q = 0; q < nb_quadrature_points; ++q) {
    std::fill_n(F, 3*3, 0.);

    for (UInt i = 0; i < spatial_dimension; ++i)
      for (UInt j = 0; j < spatial_dimension; ++j)
	F[3*i + j] = strain_val[spatial_dimension * i + j] + (i == j);

    Math::matMul<true,false>(3, 3, 3, 1., F, F, 0., C);
    Math::inv3(C, Cinv);

    Real detC = Math::det3(C);
    Real defvol = .5 * log(detC);
    Real p = lambda * defvol;
    Real traceCinv = Cinv[0] + Cinv[4] + Cinv[8];
    Real coef = mu - p;
    Real rhot = rho/sqrt(detC);
    Real tmpcele = sqrt((lambda + 2*coef)*traceCinv*traceCinv/rhot);

    cele = std::max(cele, tmpcele);
    strain_val += size_strain;
  }

  return cele;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
bool MaterialNeohookean<spatial_dimension>::setParam(const std::string & key, const std::string & value,
				  const ID & id) {
  std::stringstream sstr(value);
  if(key == "E") { sstr >> E; }
  else if(key == "nu") { sstr >> nu; }
  else if(key == "Plane_Stress") { sstr >> plane_stress; }
  else { return Material::setParam(key, value, id); }
  return true;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialNeohookean<spatial_dimension>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "MaterialNeohookean [" << std::endl;
  if(!plane_stress)
    stream << space << " + Plane strain" << std::endl;
  else
    stream << space << " + Plane stress" << std::endl;
  stream << space << " + id                      : " << id << std::endl;
  stream << space << " + name                    : " << name << std::endl;
  stream << space << " + density                 : " << rho << std::endl;
  stream << space << " + Young's modulus         : " << E << std::endl;
  stream << space << " + Poisson's ratio         : " << nu << std::endl;
  if(this->isInit()) {
    stream << space << " + First Lamé coefficient  : " << lambda << std::endl;
    stream << space << " + Second Lamé coefficient : " << mu << std::endl;
    stream << space << " + Bulk coefficient        : " << kpa << std::endl;
  }
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialNeohookean);


__END_AKANTU__
