/**
 * @file   dumper_iohelper.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Oct 26 21:52:40 2012
 *
 * @brief  implementation of DumperIOHelper
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
#include "dumper_iohelper.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
DumperIOHelper::DumperIOHelper() : count(0) { }

/* -------------------------------------------------------------------------- */
DumperIOHelper::~DumperIOHelper() {
  for (Fields::iterator it = fields.begin(); it != fields.end(); ++it) {
    delete it->second;
  }

  delete dumper;
}

/* -------------------------------------------------------------------------- */
void DumperIOHelper::dump() {
  std::stringstream filename_sstr;
  filename_sstr << filename << "_" << std::setw(4) << std::setfill('0') << count;

  try {
    dumper->dump(filename_sstr.str());
  } catch (iohelper::IOHelperException & e) {
    AKANTU_DEBUG_ERROR("I was not able to dump your data with a Dumper: " << e.what());
  }

  ++count;
}


void DumperIOHelper::registerMesh(const Mesh & mesh,
				  UInt spatial_dimension,
				  const GhostType & ghost_type,
				  const ElementKind & element_kind) {
  registerField("connectivities",
		new DumperIOHelper::ElementalField<UInt>(mesh.getConnectivities(),
							 spatial_dimension,
							 ghost_type,
							 element_kind));
  registerField("element_type",
		new DumperIOHelper::ElementTypeField(mesh,
						     spatial_dimension,
						     ghost_type,
						     element_kind));
  registerField("positions",
		new DumperIOHelper::NodalField<Real>(mesh.getNodes()));
}

/* -------------------------------------------------------------------------- */
void DumperIOHelper::registerField(const std::string & field_id,
				   Field * field) {
  Fields::iterator it = fields.find(field_id);
  if(it != fields.end())
    AKANTU_DEBUG_ERROR("The field " << field_id << " is already registered in this Dumper");

  fields[field_id] = field;
  field->registerToDumper(field_id, *dumper);
}

/* -------------------------------------------------------------------------- */
void DumperIOHelper::unRegisterField(const std::string & field_id) {
  Fields::iterator it = fields.find(field_id);
  if(it == fields.end())
    AKANTU_DEBUG_ERROR("The field " << field_id << " is not registered in this Dumper");

  delete it->second;
  fields.erase(it);
}


/* -------------------------------------------------------------------------- */
template <ElementType type>
iohelper::ElemType DumperIOHelper::getIOHelperType() {
  AKANTU_DEBUG_TO_IMPLEMENT();
  return iohelper::MAX_ELEM_TYPE;
}

template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_segment_2>() { return iohelper::LINE1; }

template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_segment_3>() { return iohelper::LINE2; }

template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_triangle_3>() { return iohelper::TRIANGLE1; }

template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_triangle_6>() { return iohelper::TRIANGLE2; }

template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_quadrangle_4>() { return iohelper::QUAD1; }

template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_quadrangle_8>() { return iohelper::QUAD2; }

template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_tetrahedron_4>() { return iohelper::TETRA1; }

template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_tetrahedron_10>() { return iohelper::TETRA2; }

template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_hexahedron_8>() { return iohelper::HEX1; }

#if defined(AKANTU_COHESIVE_ELEMENT)
template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_cohesive_2d_4>() { return iohelper::COH2D4; }

template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_cohesive_2d_6>() { return iohelper::COH2D6; }
#endif

/* -------------------------------------------------------------------------- */
iohelper::ElemType DumperIOHelper::getIOHelperType(ElementType type) {
  iohelper::ElemType ioh_type = iohelper::MAX_ELEM_TYPE;
#define GET_IOHELPER_TYPE(type)			\
  ioh_type = DumperIOHelper::getIOHelperType<type>();

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_IOHELPER_TYPE);
#undef GET_IOHELPER_TYPE
  return ioh_type;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
