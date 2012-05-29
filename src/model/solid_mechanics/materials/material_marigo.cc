/**
 * @file   material_marigo.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Chambart <marion.chambart@epfl.ch>
 * @date   Tue Jul 27 11:53:52 2010
 *
 * @brief  Specialization of the material class for the marigo material
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
#include "material_marigo.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialMarigo<spatial_dimension>::MaterialMarigo(SolidMechanicsModel & model,
						  const ID & id)  :
  Material(model, id),
  MaterialElastic<spatial_dimension>(model, id),
  MaterialDamage<spatial_dimension>(model, id),
  Yd_rand("Yd_rand",id) {
  AKANTU_DEBUG_IN();

  Yd  = 50;
  Sd  = 5000;
  Yd_randomness = 0;

  epsilon_c = std::numeric_limits<Real>::max();

  this->initInternalVector(this->Yd_rand, 1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMarigo<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialDamage<spatial_dimension>::initMaterial();

  if (std::abs(epsilon_c - std::numeric_limits<Real>::max()) < std::numeric_limits<Real>::epsilon())
    Yc = std::numeric_limits<Real>::max();
  else {
    Yc = .5 * epsilon_c * this->E * epsilon_c;
  }

  this->resizeInternalVector(this->Yd_rand);

  const Mesh & mesh = this->model->getFEM().getMesh();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension);

  for(; it != last_type; ++it) {
    UInt nb_element  = this->element_filter(*it).getSize();
    UInt nb_quad = this->model->getFEM().getNbQuadraturePoints(*it);

    Vector <Real> & Yd_rand_vec = Yd_rand(*it);

    for(UInt e = 0; e < nb_element; ++e) {
      Real rand_part = (2 * drand48()-1) * Yd_randomness * Yd;

      for(UInt q = 0; q < nb_quad; ++q)
	Yd_rand_vec(nb_quad*e+q,0) = Yd + rand_part;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMarigo<spatial_dimension>::computeStress(ElementType el_type,
						      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real F[3*3];
  Real sigma[3*3];
  Real * dam = this->damage(el_type, ghost_type).storage();
  Real * Ydq = Yd_rand(el_type, ghost_type).storage();


  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;
  memset(F, 0, 3 * 3 * sizeof(Real));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      F[3*i + j] = strain_val[spatial_dimension * i + j];

  Real Y = 0;
  computeStress(F, sigma, *dam, Y, *Ydq);
  ++dam;
  ++Ydq;

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      stress_val[spatial_dimension*i + j] = sigma[3 * i + j];

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;


  if(!this->is_non_local) this->updateDissipatedEnergy(ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
bool MaterialMarigo<spatial_dimension>::setParam(const std::string & key,
						 const std::string & value,
			       const ID & id) {
  std::stringstream sstr(value);
  if(key == "Yd") { sstr >> Yd; }
  else if(key == "Sd") { sstr >> Sd; }
  else if(key == "Yd_randomness") { sstr >> Yd_randomness; }
  else if(key == "epsilon_c") { sstr >> epsilon_c; }
  else { return MaterialDamage<spatial_dimension>::setParam(key, value, id); }
  return true;
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialMarigo<spatial_dimension>::printself(std::ostream & stream,
						  int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_marigo> [" << std::endl;
  stream << space << " + Yd         : " << Yd << std::endl;
  stream << space << " + Sd         : " << Sd << std::endl;
  stream << space << " + randomness : " << Yd_randomness << std::endl;
  if (std::abs(epsilon_c - std::numeric_limits<Real>::max()) > std::numeric_limits<Real>::epsilon()) {
    stream << space << " + epsilon_c  : " << epsilon_c << std::endl;
    stream << space << " + Yc         : " << Yc << std::endl;
  }
  MaterialDamage<spatial_dimension>::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialMarigo);


__END_AKANTU__
