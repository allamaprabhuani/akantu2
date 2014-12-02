/**
 * @file   boundary_functions.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Jun 17 14:15:57 2013
 *
 * @brief  functions for boundaries
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
// akantu
#include "aka_common.hh"
#include "solid_mechanics_model.hh"

// simtools
#include "ast_common.hh"

__BEGIN_SIMTOOLS__

using namespace akantu;

Real integrateResidual(const std::string & sub_boundary_name,
		       const SolidMechanicsModel & model,
		       UInt dir);

/// this is a fix so that all subboundaries exist on all procs
void boundaryFix(Mesh & mesh,
		 const std::vector<std::string> & sub_boundary_names);

__END_SIMTOOLS__
