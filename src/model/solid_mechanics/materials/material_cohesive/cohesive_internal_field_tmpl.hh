/**
 * @file   cohesive_internal_field_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Nov  5 22:18:00 2013
 *
 * @brief  implementation of the cohesive internal field
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

#ifndef __AKANTU_COHESIVE_INTERNAL_FIELD_TMPL_HH__
#define __AKANTU_COHESIVE_INTERNAL_FIELD_TMPL_HH__

__BEGIN_AKANTU__

template<typename T>
CohesiveInternalField<T>::CohesiveInternalField(const ID & id, Material & material) :
  InternalField<T>(id, material, material.getModel().getFEM("CohesiveFEM"),
		   dynamic_cast<MaterialCohesive &>(material).getElementFilter()) {
  this->element_kind = _ek_cohesive;

}

template<typename T>
CohesiveInternalField<T>::~CohesiveInternalField() { };

template<typename T>
void CohesiveInternalField<T>::initialize(UInt nb_component) {
  this->internalInitialize(nb_component);
}

/* -------------------------------------------------------------------------- */
template<typename T>
FacetInternalField<T>::FacetInternalField(const ID & id, Material & material) :
  InternalField<T>(id, material, material.getModel().getFEM("FacetsFEM"),
		   dynamic_cast<MaterialCohesive &>(material).getFacetFilter()) {
  this->spatial_dimension -= 1;
  this->element_kind = _ek_regular;
}

template<typename T>
FacetInternalField<T>::~FacetInternalField() { };

template<typename T>
void FacetInternalField<T>::initialize(UInt nb_component) {
  this->internalInitialize(nb_component);
}

/* -------------------------------------------------------------------------- */
template<>
inline void ParsableParamTyped< RandomInternalField<Real, FacetInternalField> >::parseParam(const ParserParameter & in_param) {
  ParsableParam::parseParam(in_param);
  RandomParameter<Real> r = in_param;
  param.setRandomDistribution(r);
}

__END_AKANTU__

#endif /* __AKANTU_COHESIVE_INTERNAL_FIELD_TMPL_HH__ */

