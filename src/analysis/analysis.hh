/**
 * @file   analysis.hh
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @date   Wed Oct  5 15:35:00 2011
 *
 * @brief  analysis object
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

#ifndef __AKANTU_ANALYSIS_HH__
#define __AKANTU_ANALYSIS_HH__

#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io_abaqus.hh"
#include "model.hh"

__BEGIN_AKANTU__


// empty policy classes
struct Abaqus {};
struct Nastran {};
// ... other possible policy classes

template <class Analysis_policy>
class Analysis;

template <>
class Analysis<Abaqus> {
  
  typedef Mesh mesh_type;
  typedef MeshIOAbaqus mesh_io_type;
  typedef Model model_type;
  
  
  mesh_type *mesh;
  model_type *model;
  
public:
  
  Analysis() : mesh(NULL), model(NULL) {}
  
  ~Analysis() {
    if (mesh) delete mesh;
    if (model) delete model;
  }
  
  /// Read Abaqus input file
  void read_file(const char* filename);
  
};




/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "analysis_inline_impl.cc"


__END_AKANTU__

#endif /* __AKANTU_ANALYSIS_HH__ */
