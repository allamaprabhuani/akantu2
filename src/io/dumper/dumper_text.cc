/**
 * @file   dumper_text.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue May 14 15:27:03 2013
 *
 * @brief  implementation of text dumper
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
#include "dumper_text.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
DumperText::DumperText(const std::string & basename, 
		       iohelper::TextDumpMode mode,
		       bool parallel) : DumperIOHelper() {
  AKANTU_DEBUG_IN();
  
  iohelper::DumperText * dumper_text = new iohelper::DumperText(mode);
  this->dumper = dumper_text;
  this->setBaseName(basename);
  
  this->setParallelContext(parallel);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DumperText::registerMesh(const Mesh & mesh,
			      UInt spatial_dimension,
			      const GhostType & ghost_type,
			      const ElementKind & element_kind) {

  registerField("positions",
		new DumperIOHelper::NodalField<Real>(mesh.getNodes()));

  // in parallel we need node type
  UInt nb_proc = StaticCommunicator::getStaticCommunicator().getNbProc();
  if (nb_proc > 1) {
    registerField("nodes_type",
		  new DumperIOHelper::NodalField<Int>(mesh.getNodesType()));
  }
}

/* -------------------------------------------------------------------------- */
void DumperText::registerFilteredMesh(const Mesh & mesh,
				      const ByElementTypeArray<UInt> & elements_filter,
				      const Array<UInt> & nodes_filter,
				      UInt spatial_dimension,
				      const GhostType & ghost_type,
				      const ElementKind & element_kind) {

  registerField("positions",
		new DumperIOHelper::NodalField<Real,
					       true>(mesh.getNodes(),
						     0,
						     0,
						     &nodes_filter));

  // in parallel we need node type
  UInt nb_proc = StaticCommunicator::getStaticCommunicator().getNbProc();
  if (nb_proc > 1) {
    registerField("nodes_type",
		  new DumperIOHelper::NodalField<Int,
						 true>(mesh.getNodesType(),
						       0,
						       0,
						       &nodes_filter));
  }
}

/* -------------------------------------------------------------------------- */
void DumperText::setBaseName(const std::string & basename) {
  AKANTU_DEBUG_IN();
  
  DumperIOHelper::setBaseName(basename);
  static_cast<iohelper::DumperText*>(this->dumper)->setDataSubDirectory(this->filename 
									+ "-DataFiles");
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DumperText::setPrecision(UInt prec) {
  AKANTU_DEBUG_IN();
  
  static_cast<iohelper::DumperText*>(this->dumper)->setPrecision(prec);

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
