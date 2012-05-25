/**
 * @file   material_vreepeerlings.cc
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 * @date   Fri Feb 17 14:00:00 2012
 *
 * @brief  Specialization of the material class for the VreePeerlings material
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
#include "material_vreepeerlings.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialVreePeerlings<spatial_dimension>::MaterialVreePeerlings(SolidMechanicsModel & model,
					     const ID & id)  :
  Material(model, id),
  MaterialElastic<spatial_dimension>(model, id),
  MaterialDamage<spatial_dimension>(model, id),
  Kapa("Kapa",id) {
  AKANTU_DEBUG_IN();

  Kapa0  = 0.0001;
  Alpha  = 0.99;
  Beta   = 300.;
  Kct    = 1.;
  Kapa0_randomness = 0.;

  this->initInternalVector(this->Kapa, 1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialVreePeerlings<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialDamage<spatial_dimension>::initMaterial();

  this->resizeInternalVector(this->Kapa);

  const Mesh & mesh = this->model->getFEM().getMesh();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension);

  for(; it != last_type; ++it) {
    Vector <Real>::iterator<Real> kapa_it  = Kapa(*it).begin();
    Vector <Real>::iterator<Real> kapa_end = Kapa(*it).end();

    for(; kapa_it != kapa_end; ++kapa_it) {
      Real rand_part = (2 * drand48()-1) * Kapa0_randomness * Kapa0;
      *kapa_it = Kapa0 + rand_part;
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialVreePeerlings<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real F[3*3];
  Real sigma[3*3];
  Real * dam = this->damage(el_type, ghost_type).storage();
  Real * Kapaq = Kapa(el_type, ghost_type).storage();


  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN;
  memset(F, 0, 3 * 3 * sizeof(Real));

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      F[3*i + j] = strain_val[spatial_dimension * i + j];

  Real Equistrain;
  computeStress(F, sigma, *dam, Equistrain, *Kapaq);
  ++dam;
  ++Kapaq;

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      stress_val[spatial_dimension*i + j] = sigma[3 * i + j];

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  if(!this->is_non_local) this->updateDissipatedEnergy(ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
bool MaterialVreePeerlings<spatial_dimension>::setParam(const std::string & key, const std::string & value,
			       const ID & id) {
  std::stringstream sstr(value);
  if(key == "Kapa0") { sstr >> Kapa0; }
  else if(key == "Alpha") { sstr >> Alpha; }
  else if(key == "Beta") { sstr >> Beta; }
  else if(key == "Kct") { sstr >> Kct; }
  else if(key == "Kapa0_randomness") { sstr >> Kapa0_randomness; }
  else { return MaterialDamage<spatial_dimension>::setParam(key, value, id); }
  return true;
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialVreePeerlings<spatial_dimension>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material<_vreepeerlings> [" << std::endl;
  stream << space << " + Kapa0            : " << Kapa0 << std::endl;
  stream << space << " + Alpha            : " << Alpha << std::endl;
  stream << space << " + Beta             : " << Beta << std::endl;
  stream << space << " + Kct              : " << Kct << std::endl;
  stream << space << " + Kapa0 randomness : " << Kapa0_randomness << std::endl;
  MaterialDamage<spatial_dimension>::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
UInt MaterialVreePeerlings<spatial_dimension>::getNbDataToPack(const Element & element,
					    SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  if(tag == _gst_smm_init_mat) {
    UInt nb_quad = this->model->getFEM().getNbQuadraturePoints(element.type);
    size += sizeof(Real) + nb_quad;
  }

  size += MaterialDamage<spatial_dimension>::getNbDataToPack(element, tag);

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
UInt MaterialVreePeerlings<spatial_dimension>::getNbDataToUnpack(const Element & element,
					      SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  if(tag == _gst_smm_init_mat) {
    UInt nb_quad = this->model->getFEM().getNbQuadraturePoints(element.type);
    size += sizeof(Real) + nb_quad;
  }

  size += MaterialDamage<spatial_dimension>::getNbDataToPack(element, tag);

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialVreePeerlings<spatial_dimension>::packData(CommunicationBuffer & buffer,
				     const Element & element,
				     SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  if(tag == _gst_smm_init_mat){
    UInt nb_quad = this->model->getFEM().getNbQuadraturePoints(element.type);
    const Vector<Real> & kapa = Kapa(element.type, _not_ghost);
    for(UInt q = 0; q < nb_quad; ++q)
      buffer << kapa(element.element * nb_quad + q);
  }

  MaterialDamage<spatial_dimension>::packData(buffer, element, tag);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialVreePeerlings<spatial_dimension>::unpackData(CommunicationBuffer & buffer,
				       const Element & element,
				       SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  if(tag == _gst_smm_init_mat) {
    UInt nb_quad = this->model->getFEM().getNbQuadraturePoints(element.type);
    Vector<Real> & kapa = Kapa(element.type, _not_ghost);
    for(UInt q = 0; q < nb_quad; ++q)
      buffer >> kapa(element.element * nb_quad + q);
  }

  MaterialDamage<spatial_dimension>::packData(buffer, element, tag);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
INSTANSIATE_MATERIAL(MaterialVreePeerlings);



__END_AKANTU__
