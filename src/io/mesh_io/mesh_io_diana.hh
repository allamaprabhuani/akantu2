/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_MESH_IO_DIANA_HH_
#define AKANTU_MESH_IO_DIANA_HH_

/* -------------------------------------------------------------------------- */
#include "mesh_io.hh"

/* -------------------------------------------------------------------------- */
#include <vector>

/* -------------------------------------------------------------------------- */

namespace akantu {

class MeshIODiana : public MeshIO {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MeshIODiana();
  ~MeshIODiana() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// read a mesh from the file
  void read(const std::string & filename, Mesh & mesh) override;

  /// write a mesh to a file
  void write(const std::string & filename, const Mesh & mesh) override;

private:
  std::string readCoordinates(std::ifstream & infile, Mesh & mesh,
                              Idx & first_node_number);

  std::string readElements(std::ifstream & infile, Mesh & mesh,
                           Idx first_node_number);

  std::string readGroups(std::ifstream & infile, Mesh & mesh,
                         Idx first_node_number);

  std::string readConnectivity(std::ifstream & infile, Mesh & mesh,
                               Idx first_node_number);

  std::string readMaterialElement(std::ifstream & infile, Mesh & mesh);

  std::string readMaterial(std::ifstream & infile,
                           const std::string & filename);

  Idx readInterval(std::stringstream & line, std::set<Idx> & interval);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  std::map<std::string, ElementType> _diana_to_akantu_element_types;
  std::map<std::string, std::string> _diana_to_akantu_mat_prop;

  /// order in witch element as to be read, akantu_node_order =
  /// _read_order[diana_node_order]
  std::map<ElementType, Int *> _read_order;

  std::map<UInt, Element> diana_element_number_to_elements;
  std::map<Element, Int> akantu_number_to_diana_number;
};

} // namespace akantu

#endif /* AKANTU_MESH_IO_DIANA_HH_ */
