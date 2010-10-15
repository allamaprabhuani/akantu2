/**
 * @file   material_elastic.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:53:52 2010
 *
 * @brief  Specialization of the material class for the elastic material
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_elastic.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialElastic::MaterialElastic(SolidMechanicsModel & model, const MaterialID & id)  :
  Material(model, id) {
  AKANTU_DEBUG_IN();

  rho = 0;
  E   = 0;
  nu  = 1./2.;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialElastic::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();

  lambda = nu * E / ((1 + nu) * (1 - 2*nu));
  mu     = E / (2 * (1 + nu));
  kpa    = lambda + 2./3. * mu;

  is_init = true;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialElastic::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real F[3*3];
  Real sigma[3*3];

  MATERIAL_QUADRATURE_POINT_LOOP_BEGIN;
  memset(F, 0, 3 * 3 * sizeof(Real));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      F[3*i + j] = strain_val[spatial_dimension * i + j];

  for (UInt i = 0; i < spatial_dimension; ++i) F[i*3 + i] -= 1;

  computeStress(F, sigma);

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      stress_val[spatial_dimension*i + j] = sigma[3 * i + j];

  MATERIAL_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialElastic::computePotentialEnergy(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if(ghost_type != _not_ghost) return;
  Real * epot = potential_energy[el_type]->values;

  MATERIAL_QUADRATURE_POINT_LOOP_BEGIN;

  computePotentialEnergy(strain_val, stress_val, epot);
  epot++;

  MATERIAL_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void MaterialElastic::setParam(const std::string & key, const std::string & value,
			       const MaterialID & id) {
  std::stringstream sstr(value);
  if(key == "rho") { sstr >> rho; }
  else if(key == "E") { sstr >> E; }
  else if(key == "nu") { sstr >> nu; }
  else { Material::setParam(key, value, id); }
}


/* -------------------------------------------------------------------------- */
void MaterialElastic::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_elastic> [" << std::endl;
  stream << space << " + id                      : " << id << std::endl;
  stream << space << " + name                    : " << name << std::endl;
  stream << space << " + density                 : " << rho << std::endl;
  stream << space << " + Young's modulus         : " << E << std::endl;
  stream << space << " + Poisson's ratio         : " << nu << std::endl;
  if(is_init) {
    stream << space << " + First Lamé coefficient  : " << lambda << std::endl;
    stream << space << " + Second Lamé coefficient : " << mu << std::endl;
    stream << space << " + Bulk coefficient        : " << kpa << std::endl;
  }
  stream << space << "]" << std::endl;
}
/* -------------------------------------------------------------------------- */

__END_AKANTU__
