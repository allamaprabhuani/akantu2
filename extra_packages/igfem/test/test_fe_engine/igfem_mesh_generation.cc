/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "mesh.hh"
#include <fstream>
/* -------------------------------------------------------------------------- */
using namespace akantu;

void generateIGFEMMesh(const ElementType type, Mesh & mesh,
                       const std::string & filename) {
  std::ifstream infile;
  infile.open(filename.c_str());
  if (!infile.good()) {
    AKANTU_ERROR("Cannot open file " << filename);
  }
  UInt current_line = 0;
  std::string line;
  UInt first_node_number = std::numeric_limits<UInt>::max();
  UInt last_node_number = 0;
  while (infile.good()) {
    std::getline(infile, line);
    current_line++;
    /// read all nodes
    if (line == "$Nodes" || line == "$NOD") {
      UInt nb_nodes;

      std::getline(infile, line);
      std::stringstream sstr(line);
      sstr >> nb_nodes;
      current_line++;

      Array<Real> & nodes = const_cast<Array<Real> &>(mesh.getNodes());
      nodes.resize(nb_nodes);

      UInt index;
      Real coord[3];
      Int spatial_dimension = nodes.getNbComponent();
      /// for each node, read the coordinates
      for (Int i = 0; i < nb_nodes; ++i) {
        UInt offset = i * spatial_dimension;

        std::getline(infile, line);
        std::stringstream sstr_node(line);
        sstr_node >> index >> coord[0] >> coord[1] >> coord[2];
        current_line++;

        first_node_number = std::min(first_node_number, index);
        last_node_number = std::max(last_node_number, index);

        /// read the coordinates
        for (Int j = 0; j < spatial_dimension; ++j)
          nodes.data()[offset + j] = coord[j];
      }
      std::getline(infile, line); /// the end of block line
    }
    /// read all elements
    if (line == "$Elements" || line == "$ELM") {
      mesh.addConnectivityType(type);
      Array<UInt> & connectivity =
          const_cast<Array<UInt> &>(mesh.getConnectivity(type));
      UInt node_per_element = connectivity.getNbComponent();
      UInt nb_elements = 0;

      std::getline(infile, line);
      std::stringstream sstr(line);
      sstr >> nb_elements;
      current_line++;

      for (Int i = 0; i < nb_elements; ++i) {
        std::getline(infile, line);
        std::stringstream sstr_elem(line);
        current_line++;
        Vector<UInt> local_connect(node_per_element);
        for (Int j = 0; j < node_per_element; ++j) {
          UInt node_index;
          sstr_elem >> node_index;

          AKANTU_DEBUG_ASSERT(node_index <= last_node_number,
                              "Node number not in range : line "
                                  << current_line);

          node_index -= first_node_number;
          local_connect(j) = node_index;
        }
        connectivity.push_back(local_connect);
      }
      std::getline(infile, line); /// the end of block line
    }
  }
}
