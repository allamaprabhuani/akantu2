/**
 * @file   material_elastic.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:53:52 2010
 *
 * @brief  Specialization of the material class for the elastic material
 *
 * @section LICENSE
 *
 * <insert license here>
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
void MaterialElastic::constitutiveLaw(ElementType el_type, bool local) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension    = model.getSpatialDimension();
  UInt nb_quadrature_points = FEM::getNbQuadraturePoints(el_type);
  UInt size_strain          = spatial_dimension * spatial_dimension;

  UInt nb_element;
  Real * strain_val;
  Real * stress_val;

#ifdef AKANTU_USE_MPI
  if(local) {
#endif //AKANTU_USE_MPI
    nb_element   = element_filter[el_type]->getSize();
    strain_val = strain[el_type]->values;
    stress_val = stress[el_type]->values;
#ifdef AKANTU_USE_MPI
  } else {
    nb_element = element_filter[el_type]->getSize();
    strain_val = strain[el_type]->values;
    stress_val = stress[el_type]->values;
    potential_energy_flag = false;
  }
#endif //AKANTU_USE_MPI

  Real * epot = NULL;
  if (potential_energy_flag) epot = potential_energy[el_type]->values;

  Real F[3*3];
  Real sigma[3*3];

  /// for each element
  for (UInt el = 0; el < nb_element; ++el) {
    /// for each quadrature points
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      memset(F, 0, 3 * 3 * sizeof(Real));
      for (UInt i = 0; i < spatial_dimension; ++i)
	for (UInt j = 0; j < spatial_dimension; ++j)
	  F[3*i + j] = strain_val[spatial_dimension * i + j];

      for (UInt i = 0; i < spatial_dimension; ++i) F[i*3 + i] -= 1;

      constitutiveLaw(F, sigma, epot);

      for (UInt i = 0; i < spatial_dimension; ++i)
	for (UInt j = 0; j < spatial_dimension; ++j)
	  stress_val[spatial_dimension*i + j] = sigma[3 * i + j];

      strain_val += size_strain;
      stress_val += size_strain;
      if (potential_energy_flag) epot += nb_quadrature_points;
    }
  }
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
