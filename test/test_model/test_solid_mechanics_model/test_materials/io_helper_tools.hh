/**
 * @file   io_helper_tools.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Sep 30 11:15:18 2011
 *
 * @brief  
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
#include "aka_common.hh"
#include "solid_mechanics_model.hh"
#include "io_helper.hh"

extern const akantu::UInt spatial_dimension;

/* ------------------------------------------------------------------------ */
iohelper::ElemType getIOHelperType(akantu::ElementType type);
/* -------------------------------------------------------------------------- */
void paraviewInit(iohelper::Dumper & dumper,
		  const akantu::SolidMechanicsModel & model,
		  const akantu::ElementType & type,
		  const std::string & filename);

void paraviewDump(iohelper::Dumper & dumper);

// void checkpointInit(iohelper::Dumper & dumper,
// 		    const akantu::SolidMechanicsModel & model,
// 		    const akantu::ElementType & type,
// 		    const std::string & filename);

// void checkpoint(iohelper::Dumper & dumper,
// 		const akantu::SolidMechanicsModel & model);

// void restart(const akantu::SolidMechanicsModel & model,
// 	     const akantu::ElementType & type,
// 	     const std::string & filename);
