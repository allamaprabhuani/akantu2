/**
 * @file   material_viscoelastic.hh
 * @author Vlad Yastrebov <vladislav.yastrebov@epfl.ch> 
 * @date   Thu Feb 7 2012 
 *
 * @brief  Material Visco-elastic
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
#include "material_viscoelastic.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialViscoElastic::MaterialViscoElastic(Model & model, const ID & id)  :
  Material(model, id), MaterialElastic(model, id),
  history_integral("history_integral", id),
  stress_dev("stress_dev", id) {
  AKANTU_DEBUG_IN();

//  memset(stress_deviator, 0, 3 * 3 * sizeof(Real));
//  memset(h, 0, 3 * 3 * sizeof(Real));

  eta = 1.;
  Ev  = 1.;
//  for (UInt i = 0; i < spatial_dimension; ++i)
//    for (UInt j = 0; j < spatial_dimension; ++j)
//      h[i*3 + j] = 0.;
  UInt spatial_dimension = this->model->getSpatialDimension();
  UInt stress_size = spatial_dimension * spatial_dimension;

  initInternalVector(this->stress_dev, stress_size);
  initInternalVector(this->history_integral, stress_size);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialViscoElastic::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialElastic::initMaterial();

  resizeInternalVector(this->stress_dev);
  resizeInternalVector(this->history_integral);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialViscoElastic::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real Theta, tau = eta/Ev, K = E / ( 3 * (1 - 2*nu) );
  Real e[3*3], s[3*3];
  Real sigma[3*3];

  Vector<Real> & stress_d  = stress_dev(el_type, ghost_type);
  Vector<Real> & history_i = history_integral(el_type, ghost_type);

  Real * stress_d_val  = stress_d.storage();
  Real * history_i_val = history_i.storage();
  Real dt = model->getTimeStep();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;
  memset(e, 0, 3 * 3 * sizeof(Real));
  memset(s, 0, 3 * 3 * sizeof(Real));
  memset(sigma, 0, 3 * 3 * sizeof(Real));

  /// Compute the first invariant of strain
  Theta = 0;
  for (UInt i = 0; i < spatial_dimension; ++i)
    Theta += strain_val[spatial_dimension * i + i];

  /// Compute the deviator of strain
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      if (i == j)
        e[3*i + j] = strain_val[spatial_dimension * i + j] - 1/double(spatial_dimension)*Theta;
      else
        e[3*i + j] = strain_val[spatial_dimension * i + j];

  /// Compute the deviator of stress
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
       s[3*i + j] = E / (1 + nu) * e[3*i + j];

  /// Update the internal variable
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      history_i_val[spatial_dimension * i + j]  = exp( -dt/tau ) * history_i_val[spatial_dimension * i + j] + exp( -0.5*dt/tau ) * ( s[3*i + j] - stress_d_val[spatial_dimension*i + j]);

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)  
      if (i == j)
        sigma[3*i + j] = 2./3. * K * Theta + ( E * s[3*i + j] + Ev * history_i_val[spatial_dimension*i + j] ) / ( Ev + E );
      else
        sigma[3*i + j] = ( E * s[3*i + j] + Ev * history_i_val[spatial_dimension*i + j] ) / ( Ev + E );
  
  /// Save the deviator of stress
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
       stress_d_val[spatial_dimension*i + j] = s[3*i + j];

//  computeStress(F, sigma);

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      stress_val[spatial_dimension*i + j] = sigma[3 * i + j];

  stress_d_val += spatial_dimension * spatial_dimension;
  history_i_val += spatial_dimension * spatial_dimension;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

/*
  UInt spatial_dimension = model->getSpatialDimension();

  Vector<UInt> & elem_filter = element_filter  (el_type, ghost_type);
  Vector<Real> & stress_visc = stress_viscosity(el_type, ghost_type);
  Vector<Real> & stress_el   = stress_elastic  (el_type, ghost_type);
  //mu       = E / (2 * (1 + nu));

  MaterialElastic::computeStress(el_type, ghost_type);

  Vector<Real> & velocity = model->getVelocity();
  UInt nb_elem = element_filter(el_type, ghost_type).getSize();
  UInt nb_quad_points_per_elem =
    model->getFEM().getNbQuadraturePoints(el_type, ghost_type);
  UInt nb_quad_points = nb_quad_points_per_elem * nb_elem;

  Vector<Real> strain_rate(nb_quad_points,
			   spatial_dimension * spatial_dimension,
			   "strain_rate");

  model->getFEM().gradientOnQuadraturePoints(velocity, strain_rate,
   					     spatial_dimension,
   					     el_type, ghost_type, &elem_filter);

  Real F[3*3];
  Real sigma[3*3];
  Real * strain_rate_val = strain_rate.storage();
  Real * stress_visc_val = stress_visc.storage();
  Real * stress_el_val   = stress_el.storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;
  memset(F, 0, 3 * 3 * sizeof(Real));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      F[3*i + j] = strain_rate_val[spatial_dimension * i + j];

  MaterialElastic::computeStress(F, sigma);

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j) {
//      stress_visc_val[spatial_dimension * i + j]  = alpha * sigma[3 * i + j];
      stress_el_val  [spatial_dimension * i + j]  = stress_val[spatial_dimension*i + j];
//      stress_val     [spatial_dimension * i + j] += stress_visc_val[spatial_dimension*i + j];
    }

  strain_rate_val += spatial_dimension * spatial_dimension;
//  stress_visc_val += spatial_dimension * spatial_dimension;
  stress_el_val   += spatial_dimension * spatial_dimension;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
*/
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialViscoElastic::computePotentialEnergy(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();
/*
  if(ghost_type != _not_ghost) return;

  Vector<Real> & stress_el = stress_elastic(el_type, ghost_type);
  Real * stress_el_val = stress_el.storage();

  Real * epot = potential_energy(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;

  MaterialElastic::computePotentialEnergy(strain_val, stress_el_val, epot);
  epot++;
  stress_el_val += spatial_dimension*spatial_dimension;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
*/
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool MaterialViscoElastic::setParam(const std::string & key, const std::string & value,
				      const ID & id) {
  std::stringstream sstr(value);
  if(key == "eta") { sstr >> eta; }
  else if(key == "Ev") { sstr >> Ev; }
  else { return MaterialElastic::setParam(key, value, id); }
  return true;
}


/* -------------------------------------------------------------------------- */
void MaterialViscoElastic::printself(std::ostream & stream, int indent) const {
/*
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "MaterialViscoElastic [" << std::endl;
  MaterialElastic::printself(stream, indent + 1);
  stream << space << " + artifical viscous ratio : " << alpha << std::endl;
  stream << space << "]" << std::endl;
*/
}
/* -------------------------------------------------------------------------- */


__END_AKANTU__
