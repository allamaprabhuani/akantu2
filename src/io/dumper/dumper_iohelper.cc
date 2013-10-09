/**
 * @file   dumper_iohelper.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author David Kammer <david.kammer@epfl.ch>
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
#include "sub_boundary.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
DumperIOHelper::DumperIOHelper() : count(0), time_activated(false) {}

/* -------------------------------------------------------------------------- */
DumperIOHelper::~DumperIOHelper() {
  for (Fields::iterator it = fields.begin(); it != fields.end(); ++it) {
    delete it->second;
  }

  delete dumper;
}

/* -------------------------------------------------------------------------- */
void DumperIOHelper::setParallelContext(bool is_parallel) {
  UInt whoami  = StaticCommunicator::getStaticCommunicator().whoAmI();
  UInt nb_proc = StaticCommunicator::getStaticCommunicator().getNbProc();

  if(is_parallel)
    dumper->setParallelContext(whoami, nb_proc);
  else
    dumper->setParallelContext(0, 1);
}

/* -------------------------------------------------------------------------- */
void DumperIOHelper::setDirectory(const std::string & directory) {
  this->directory = directory;
  dumper->setPrefix(directory);
}

/* -------------------------------------------------------------------------- */
void DumperIOHelper::setBaseName(const std::string & basename) {
  filename = basename;
}

/* -------------------------------------------------------------------------- */
void DumperIOHelper::setTimeStep(Real time_step) {
  if(!time_activated)
    this->dumper->activateTimeDescFiles(time_step);
  else
    this->dumper->setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
void DumperIOHelper::dump() {
  try {
    dumper->dump(filename, count);
  } catch (iohelper::IOHelperException & e) {
    AKANTU_DEBUG_ERROR("I was not able to dump your data with a Dumper: " << e.what());
  }

  ++count;
}

/* -------------------------------------------------------------------------- */
void DumperIOHelper::dump(UInt step) {
  this->count = step;
  this->dump();
}

/* -------------------------------------------------------------------------- */
void DumperIOHelper::dump(Real current_time, UInt step) {
  this->dumper->setCurrentTime(current_time);
  this->dump(step);
}

/* -------------------------------------------------------------------------- */
void DumperIOHelper::registerMesh(const Mesh & mesh,
				  UInt spatial_dimension,
				  const GhostType & ghost_type,
				  const ElementKind & element_kind) {
#if defined(AKANTU_COHESIVE_ELEMENT)
  if (element_kind == _ek_cohesive) {
    registerField("connectivities",
		  new DumperIOHelper::CohesiveConnectivityField(mesh,
								spatial_dimension,
								ghost_type));
  } else
#endif
    registerField("connectivities",
		  new DumperIOHelper::ElementalField<UInt>(mesh.getConnectivities(),
							   spatial_dimension,
							   ghost_type,
							   element_kind));

  registerField("element_type",
		new DumperIOHelper::ElementTypeField<>(mesh,
						     spatial_dimension,
						     ghost_type,
						     element_kind));
  registerField("positions",
		new DumperIOHelper::NodalField<Real>(mesh.getNodes()));
}

/* -------------------------------------------------------------------------- */
void DumperIOHelper::registerFilteredMesh(const Mesh & mesh,
					  const ByElementTypeArray<UInt> & elements_filter,
					  const Array<UInt> & nodes_filter,
					  UInt spatial_dimension,
					  const GhostType & ghost_type,
					  const ElementKind & element_kind) {
  registerField("connectivities",
		new DumperIOHelper::FilteredConnectivityField(mesh.getConnectivities(),
							      elements_filter,
							      nodes_filter,
							      spatial_dimension,
							      ghost_type,
							      element_kind));
  registerField("element_type",
		new DumperIOHelper::ElementTypeField<true>(mesh,
							   spatial_dimension,
							   ghost_type,
							   element_kind,
							   &elements_filter));

  registerField("positions",
		new DumperIOHelper::NodalField<Real,
					       true>(mesh.getNodes(),
						     0,
						     0,
						     &nodes_filter));
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
void DumperIOHelper::registerVariable(const std::string & variable_id,
				      VariableBase * variable) {
  Variables::iterator it = variables.find(variable_id);
  if(it != variables.end())
    AKANTU_DEBUG_ERROR("The Variable " << variable_id << " is already registered in this Dumper");

  variables[variable_id] = variable;
  variable->registerToDumper(variable_id, *dumper);
}

/* -------------------------------------------------------------------------- */
void DumperIOHelper::unRegisterVariable(const std::string & variable_id) {
  Variables::iterator it = variables.find(variable_id);
  if(it == variables.end())
    AKANTU_DEBUG_ERROR("The variable " << variable_id << " is not registered in this Dumper");

  delete it->second;
  variables.erase(it);
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

template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_pentahedron_6>() { return iohelper::PRISM1; }

#if defined(AKANTU_COHESIVE_ELEMENT)
template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_cohesive_2d_4>() { return iohelper::COH2D4; }

template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_cohesive_2d_6>() { return iohelper::COH2D6; }

template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_cohesive_3d_6>() { return iohelper::COH3D6; }

template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_cohesive_3d_12>() { return iohelper::COH3D12; }
#endif

#if defined(AKANTU_STRUCTURAL_MECHANICS)
template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_bernoulli_beam_2>() { return iohelper::LINE1; }

template <>
iohelper::ElemType DumperIOHelper::getIOHelperType<_bernoulli_beam_3>() { return iohelper::LINE2; }
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
