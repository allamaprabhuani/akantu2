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

/* -----------------------------------------------------------------------------
   Version (Legacy) 1.0

   $NOD
   number-of-nodes
   node-number x-coord y-coord z-coord
   ...
   $ENDNOD
   $ELM
   number-of-elements
   elm-number elm-type reg-phys reg-elem number-of-nodes node-number-list
   ...
   $ENDELM
   -----------------------------------------------------------------------------
   Version 2.1

   $MeshFormat
   version-number file-type data-size
   $EndMeshFormat
   $Nodes
   number-of-nodes
   node-number x-coord y-coord z-coord
   ...
   $EndNodes
   $Elements
   number-of-elements
   elm-number elm-type number-of-tags < tag > ... node-number-list
   ...
   $EndElements
   $PhysicalNames
   number-of-names
   physical-dimension physical-number "physical-name"
   ...
   $EndPhysicalNames
   $NodeData
   number-of-string-tags
   < "string-tag" >
   ...
   number-of-real-tags
   < real-tag >
   ...
   number-of-integer-tags
   < integer-tag >
   ...
   node-number value ...
   ...
   $EndNodeData
   $ElementData
   number-of-string-tags
   < "string-tag" >
   ...
   number-of-real-tags
   < real-tag >
   ...
   number-of-integer-tags
   < integer-tag >
   ...
   elm-number value ...
   ...
   $EndElementData
   $ElementNodeData
   number-of-string-tags
   < "string-tag" >
   ...
   number-of-real-tags
   < real-tag >
   ...
   number-of-integer-tags
   < integer-tag >
   ...
   elm-number number-of-nodes-per-element value ...
   ...
   $ElementEndNodeData

   -----------------------------------------------------------------------------
   elem-type

   1:  2-node line.
   2:  3-node triangle.
   3:  4-node quadrangle.
   4:  4-node tetrahedron.
   5:  8-node hexahedron.
   6:  6-node prism.
   7:  5-node pyramid.
   8:  3-node second order line
   9:  6-node second order triangle
   10: 9-node second order quadrangle
   11: 10-node second order tetrahedron
   12: 27-node second order hexahedron
   13: 18-node second order prism
   14: 14-node second order pyramid
   15: 1-node point.
   16: 8-node second order quadrangle
   17: 20-node second order hexahedron
   18: 15-node second order prism
   19: 13-node second order pyramid
   20: 9-node third order incomplete triangle
   21: 10-node third order triangle
   22: 12-node fourth order incomplete triangle
   23: 15-node fourth order triangle
   24: 15-node fifth order incomplete triangle
   25: 21-node fifth order complete triangle
   26: 4-node third order edge
   27: 5-node fourth order edge
   28: 6-node fifth order edge
   29: 20-node third order tetrahedron
   30: 35-node fourth order tetrahedron
   31: 56-node fifth order tetrahedron
   -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#include "element_group.hh"
#include "mesh_io.hh"
#include "mesh_utils.hh"
#include "node_group.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
/*   Methods Implentations                                                    */
/* -------------------------------------------------------------------------- */
MeshIOMSH::MeshIOMSH() {
  canReadSurface = true;
  canReadExtendedData = true;

  _msh_nodes_per_elem[_msh_not_defined] = 0;
  _msh_nodes_per_elem[_msh_segment_2] = 2;
  _msh_nodes_per_elem[_msh_triangle_3] = 3;
  _msh_nodes_per_elem[_msh_quadrangle_4] = 4;
  _msh_nodes_per_elem[_msh_tetrahedron_4] = 4;
  _msh_nodes_per_elem[_msh_hexahedron_8] = 8;
  _msh_nodes_per_elem[_msh_prism_1] = 6;
  _msh_nodes_per_elem[_msh_pyramid_1] = 1;
  _msh_nodes_per_elem[_msh_segment_3] = 3;
  _msh_nodes_per_elem[_msh_triangle_6] = 6;
  _msh_nodes_per_elem[_msh_quadrangle_9] = 9;
  _msh_nodes_per_elem[_msh_tetrahedron_10] = 10;
  _msh_nodes_per_elem[_msh_hexahedron_27] = 27;
  _msh_nodes_per_elem[_msh_hexahedron_20] = 20;
  _msh_nodes_per_elem[_msh_prism_18] = 18;
  _msh_nodes_per_elem[_msh_prism_15] = 15;
  _msh_nodes_per_elem[_msh_pyramid_14] = 14;
  _msh_nodes_per_elem[_msh_point] = 1;
  _msh_nodes_per_elem[_msh_quadrangle_8] = 8;

  _msh_to_akantu_element_types[_msh_not_defined] = _not_defined;
  _msh_to_akantu_element_types[_msh_segment_2] = _segment_2;
  _msh_to_akantu_element_types[_msh_triangle_3] = _triangle_3;
  _msh_to_akantu_element_types[_msh_quadrangle_4] = _quadrangle_4;
  _msh_to_akantu_element_types[_msh_tetrahedron_4] = _tetrahedron_4;
  _msh_to_akantu_element_types[_msh_hexahedron_8] = _hexahedron_8;
  _msh_to_akantu_element_types[_msh_prism_1] = _pentahedron_6;
  _msh_to_akantu_element_types[_msh_pyramid_1] = _not_defined;
  _msh_to_akantu_element_types[_msh_segment_3] = _segment_3;
  _msh_to_akantu_element_types[_msh_triangle_6] = _triangle_6;
  _msh_to_akantu_element_types[_msh_quadrangle_9] = _not_defined;
  _msh_to_akantu_element_types[_msh_tetrahedron_10] = _tetrahedron_10;
  _msh_to_akantu_element_types[_msh_hexahedron_27] = _not_defined;
  _msh_to_akantu_element_types[_msh_hexahedron_20] = _hexahedron_20;
  _msh_to_akantu_element_types[_msh_prism_18] = _not_defined;
  _msh_to_akantu_element_types[_msh_prism_15] = _pentahedron_15;
  _msh_to_akantu_element_types[_msh_pyramid_14] = _not_defined;
  _msh_to_akantu_element_types[_msh_point] = _point_1;
  _msh_to_akantu_element_types[_msh_quadrangle_8] = _quadrangle_8;

  _akantu_to_msh_element_types[_not_defined] = _msh_not_defined;
  _akantu_to_msh_element_types[_segment_2] = _msh_segment_2;
  _akantu_to_msh_element_types[_segment_3] = _msh_segment_3;
  _akantu_to_msh_element_types[_triangle_3] = _msh_triangle_3;
  _akantu_to_msh_element_types[_triangle_6] = _msh_triangle_6;
  _akantu_to_msh_element_types[_tetrahedron_4] = _msh_tetrahedron_4;
  _akantu_to_msh_element_types[_tetrahedron_10] = _msh_tetrahedron_10;
  _akantu_to_msh_element_types[_quadrangle_4] = _msh_quadrangle_4;
  _akantu_to_msh_element_types[_quadrangle_8] = _msh_quadrangle_8;
  _akantu_to_msh_element_types[_hexahedron_8] = _msh_hexahedron_8;
  _akantu_to_msh_element_types[_hexahedron_20] = _msh_hexahedron_20;
  _akantu_to_msh_element_types[_pentahedron_6] = _msh_prism_1;
  _akantu_to_msh_element_types[_pentahedron_15] = _msh_prism_15;
  _akantu_to_msh_element_types[_point_1] = _msh_point;

#if defined(AKANTU_STRUCTURAL_MECHANICS)
  _akantu_to_msh_element_types[_bernoulli_beam_2] = _msh_segment_2;
  _akantu_to_msh_element_types[_bernoulli_beam_3] = _msh_segment_2;
  _akantu_to_msh_element_types[_discrete_kirchhoff_triangle_18] =
      _msh_triangle_3;
#endif

  std::map<ElementType, MSHElementType>::iterator it;
  for (it = _akantu_to_msh_element_types.begin();
       it != _akantu_to_msh_element_types.end(); ++it) {
    Int nb_nodes = _msh_nodes_per_elem[it->second];

    std::vector<Idx> tmp(nb_nodes);
    for (Int i = 0; i < nb_nodes; ++i) {
      tmp[i] = i;
    }

    switch (it->first) {
    case _tetrahedron_10:
      tmp[8] = 9;
      tmp[9] = 8;
      break;
    case _pentahedron_6:
      tmp[0] = 2;
      tmp[1] = 0;
      tmp[2] = 1;
      tmp[3] = 5;
      tmp[4] = 3;
      tmp[5] = 4;
      break;
    case _pentahedron_15:
      tmp[0] = 2;
      tmp[1] = 0;
      tmp[2] = 1;
      tmp[3] = 5;
      tmp[4] = 3;
      tmp[5] = 4;
      tmp[6] = 8;
      tmp[8] = 11;
      tmp[9] = 6;
      tmp[10] = 9;
      tmp[11] = 10;
      tmp[12] = 14;
      tmp[14] = 12;
      break;
    case _hexahedron_20:
      tmp[9] = 11;
      tmp[10] = 12;
      tmp[11] = 9;
      tmp[12] = 13;
      tmp[13] = 10;
      tmp[17] = 19;
      tmp[18] = 17;
      tmp[19] = 18;
      break;
    default:
      // nothing to change
      break;
    }
    _read_order[it->first] = tmp;
  }
}

/* -------------------------------------------------------------------------- */
MeshIOMSH::~MeshIOMSH() = default;

/* -------------------------------------------------------------------------- */
namespace {
  struct File {
    std::string filename;
    std::ifstream infile;
    std::string line;
    std::size_t current_line{0};

    std::size_t first_node_number{std::numeric_limits<std::size_t>::max()},
        last_node_number{0};
    bool has_physical_names{false};

    std::unordered_map<size_t, size_t> node_tags;
    std::unordered_map<size_t, Element> element_tags;
    int version{0};
    int size_of_size_t{0};
    Mesh & mesh;
    MeshAccessor mesh_accessor;

    std::multimap<std::pair<int, int>, int> entity_tag_to_physical_tags;

    File(const std::string & filename, Mesh & mesh)
        : filename(filename), mesh(mesh), mesh_accessor(mesh) {
      infile.open(filename.c_str());
      if (not infile.good()) {
        AKANTU_EXCEPTION("Cannot open file " << filename);
      }
    }

    ~File() { infile.close(); }

    auto good() { return infile.good(); }

    std::stringstream get_line() {
      std::string tmp_str;
      if (infile.eof()) {
        AKANTU_EXCEPTION("Reached the end of the file " << filename);
      }
      std::getline(infile, tmp_str);
      line = trim(tmp_str);
      ++current_line;

      return std::stringstream(line);
    }

    template <typename... Ts> void read_line(Ts &&... ts) {
      auto && sstr = get_line();
      (void)std::initializer_list<int>{
          (sstr >> std::forward<decltype(ts)>(ts), 0)...};
    }
  };
} // namespace

/* -------------------------------------------------------------------------- */
template <typename File, typename Readers>
void MeshIOMSH::populateReaders2(File & file, Readers & readers) {
  readers["$NOD"] = readers["$Nodes"] = [&](const std::string & /*unused*/) {
    UInt nb_nodes;
    file.read_line(nb_nodes);

    Array<Real> & nodes = file.mesh_accessor.getNodes();
    nodes.resize(nb_nodes);
    file.mesh_accessor.setNbGlobalNodes(nb_nodes);

    size_t index;
    Vector<double> coord(3);

    /// for each node, read the coordinates
    for (auto && data : enumerate(make_view(nodes, nodes.getNbComponent()))) {
      file.read_line(index, coord(0), coord(1), coord(2));

      if (index > std::numeric_limits<UInt>::max()) {
        AKANTU_EXCEPTION(
            "There are more nodes in this files than the index type of akantu "
            "can handle, consider recompiling with a bigger index type");
      }

      file.first_node_number = std::min(file.first_node_number, index);
      file.last_node_number = std::max(file.last_node_number, index);

      for (auto && coord_data : zip(std::get<1>(data), coord)) {
        std::get<0>(coord_data) = std::get<1>(coord_data);
      }
      file.node_tags[index] = std::get<0>(data);
    }
  };

  readers["$ELM"] = readers["$Elements"] = [&](const std::string & /*unused*/) {
    Int nb_elements;
    file.read_line(nb_elements);

    Int index;
    UInt msh_type;
    ElementType akantu_type;

    for (Int i = 0; i < nb_elements; ++i) {
      auto && sstr_elem = file.get_line();

      sstr_elem >> index;
      sstr_elem >> msh_type;

      /// get the connectivity vector depending on the element type
      akantu_type =
          this->_msh_to_akantu_element_types[MSHElementType(msh_type)];

      if (akantu_type == _not_defined) {
        AKANTU_DEBUG_WARNING("Unsuported element kind "
                             << msh_type << " at line " << file.current_line);
        continue;
      }

      Element elem{akantu_type, 0, _not_ghost};

      auto & connectivity = file.mesh_accessor.getConnectivity(akantu_type);
      auto node_per_element = connectivity.getNbComponent();
      auto & read_order = this->_read_order[akantu_type];

      /// read tags informations
      if (file.version < 2000) {
        Int tag0;
        Int tag1;
        Int nb_nodes; // reg-phys, reg-elem, number-of-nodes
        sstr_elem >> tag0 >> tag1 >> nb_nodes;

        auto & data0 =
            file.mesh_accessor.template getData<Int>("tag_0", akantu_type);
        data0.push_back(tag0);

        auto & data1 =
            file.mesh_accessor.template getData<Int>("tag_1", akantu_type);
        data1.push_back(tag1);
      } else if (file.version < 4000) {
        Int nb_tags;
        sstr_elem >> nb_tags;
        for (Int j = 0; j < nb_tags; ++j) {
          Int tag;
          sstr_elem >> tag;

          auto & data = file.mesh_accessor.template getData<Int>(
              "tag_" + std::to_string(j), akantu_type);
          data.push_back(tag);
        }
      }

      Vector<Idx> local_connect(node_per_element);
      for (Int j = 0; j < node_per_element; ++j) {
        Int node_index;
        sstr_elem >> node_index;

        AKANTU_DEBUG_ASSERT(node_index <= Int(file.last_node_number),
                            "Node number not in range : line "
                                << file.current_line);

        local_connect(read_order[j]) = file.node_tags[node_index];
      }

      connectivity.push_back(local_connect);
      elem.element = connectivity.size() - 1;
      file.element_tags[index] = elem;
    }
  };

  readers["$Periodic"] = [&](const std::string &) {
    Int nb_periodic_entities;
    file.read_line(nb_periodic_entities);

    file.mesh_accessor.getNodesFlags().resize(file.mesh.getNbNodes(),
                                              NodeFlag::_normal);

    for (Int p = 0; p < nb_periodic_entities; ++p) {
      // dimension slave-tag master-tag
      Int dimension;
      file.read_line(dimension);

      // transformation
      file.get_line();

      // nb nodes
      Int nb_nodes;
      file.read_line(nb_nodes);

      for (Int n = 0; n < nb_nodes; ++n) {
        // slave master
        auto && sstr = file.get_line();

        // The info in the mesh seem inconsistent so they are ignored for now.
        continue;

        if (dimension == file.mesh.getSpatialDimension() - 1) {
          Idx slave, master;

          sstr >> slave;
          sstr >> master;
          file.mesh_accessor.addPeriodicSlave(file.node_tags[slave],
                                              file.node_tags[master]);
        }
      }
    }

    // mesh_accessor.markMeshPeriodic();
  };
}

/* -------------------------------------------------------------------------- */
template <typename File, typename Readers>
void MeshIOMSH::populateReaders4(File & file, Readers & readers) {
  static std::map<int, std::string> entity_type{
      {0, "points"},
      {1, "curve"},
      {2, "surface"},
      {3, "volume"},
  };

  readers["$Entities"] = [&](const std::string & /*unused*/) {
    size_t num_entity[4];
    file.read_line(num_entity[0], num_entity[1], num_entity[2], num_entity[3]);

    for (auto entity_dim : arange(4)) {
      for (auto _ [[gnu::unused]] : arange(num_entity[entity_dim])) {
        auto && sstr = file.get_line();

        int tag;
        double min_x;
        double min_y;
        double min_z;
        double max_x;
        double max_y;
        double max_z;
        size_t num_physical_tags;
        sstr >> tag >> min_x >> min_y >> min_z;

        if (entity_dim > 0 or file.version < 4001) {
          sstr >> max_x >> max_y >> max_z;
        }

        sstr >> num_physical_tags;

        for (auto _ [[gnu::unused]] : arange(num_physical_tags)) {
          int phys_tag;
          sstr >> phys_tag;

          std::string physical_name;
          if (this->physical_names.find(phys_tag) ==
              this->physical_names.end()) {
            physical_name = "msh_block_" + std::to_string(phys_tag);
          } else {
            physical_name = this->physical_names[phys_tag];
          }

          if (not file.mesh.elementGroupExists(physical_name)) {
            file.mesh.createElementGroup(physical_name, entity_dim);
          } else {
            file.mesh.getElementGroup(physical_name).addDimension(entity_dim);
          }
          file.entity_tag_to_physical_tags.insert(
              std::make_pair(std::make_pair(tag, entity_dim), phys_tag));
        }
      }
    }
  };

  readers["$Nodes"] = [&](const std::string & /*unused*/) {
    size_t num_blocks;
    size_t num_nodes;
    if (file.version >= 4001) {
      file.read_line(num_blocks, num_nodes, file.first_node_number,
                     file.last_node_number);
    } else {
      file.read_line(num_blocks, num_nodes);
    }
    auto & nodes = file.mesh_accessor.getNodes();
    nodes.reserve(num_nodes);
    file.mesh_accessor.setNbGlobalNodes(num_nodes);

    if (num_nodes > std::numeric_limits<UInt>::max()) {
      AKANTU_EXCEPTION(
          "There are more nodes in this files than the index type of akantu "
          "can handle, consider recompiling with a bigger index type");
    }

    size_t node_id{0};

    auto dim = nodes.getNbComponent();

    for (auto block [[gnu::unused]] : arange(num_blocks)) {
      int entity_dim;
      int entity_tag;
      int parametric;
      size_t num_nodes_in_block;
      Vector<double> pos(3);

      if (file.version >= 4001) {
        file.read_line(entity_dim, entity_tag, parametric, num_nodes_in_block);
        if (parametric) {
          AKANTU_EXCEPTION(
              "Akantu does not support parametric nodes in msh files");
        }
        for (auto _ [[gnu::unused]] : arange(num_nodes_in_block)) {
          size_t tag;
          file.read_line(tag);
          file.node_tags[tag] = node_id;
          ++node_id;
        }

        for (auto _ [[gnu::unused]] : arange(num_nodes_in_block)) {
          file.read_line(pos(_x), pos(_y), pos(_z));
          nodes.push_back(pos.block(0, 0, dim, 1));
        }
      } else {
        file.read_line(entity_tag, entity_dim, parametric, num_nodes_in_block);

        for (auto _ [[gnu::unused]] : arange(num_nodes_in_block)) {
          size_t tag;
          file.read_line(tag, pos(_x), pos(_y), pos(_z));

          file.first_node_number = std::min(file.first_node_number, tag);
          file.last_node_number = std::max(file.last_node_number, tag);

          nodes.push_back(pos.block(0, 0, dim, 1));

          file.node_tags[tag] = node_id;
          ++node_id;
        }
      }
    }
  };

  readers["$Elements"] = [&](const std::string & /*unused*/) {
    size_t num_blocks;
    size_t num_elements;
    file.read_line(num_blocks, num_elements);

    for (auto block [[gnu::unused]] : arange(num_blocks)) {
      int entity_dim;
      int entity_tag;
      int element_type;
      size_t num_elements_in_block;

      if (file.version >= 4001) {
        file.read_line(entity_dim, entity_tag, element_type,
                       num_elements_in_block);
      } else {
        file.read_line(entity_tag, entity_dim, element_type,
                       num_elements_in_block);
      }

      /// get the connectivity vector depending on the element type
      auto && akantu_type =
          this->_msh_to_akantu_element_types[(MSHElementType)element_type];

      if (akantu_type == _not_defined) {
        AKANTU_DEBUG_WARNING("Unsuported element kind " << element_type
                                                        << " at line "
                                                        << file.current_line);
        continue;
      }

      Element elem{akantu_type, 0, _not_ghost};

      auto & connectivity = file.mesh_accessor.getConnectivity(akantu_type);
      Vector<Idx> local_connect(connectivity.getNbComponent());
      auto && read_order = this->_read_order[akantu_type];

      auto & data0 =
          file.mesh_accessor.template getData<Idx>("tag_0", akantu_type);
      data0.resize(data0.size() + num_elements_in_block, 0);

      auto range = file.entity_tag_to_physical_tags.equal_range(
          std::make_pair(entity_tag, entity_dim));

      auto & physical_data = file.mesh_accessor.template getData<std::string>(
          "physical_names", akantu_type);
      physical_data.resize(physical_data.size() + num_elements_in_block, "");

      for (auto _ [[gnu::unused]] : arange(num_elements_in_block)) {
        auto && sstr_elem = file.get_line();
        std::size_t elem_tag;
        sstr_elem >> elem_tag;
        for (auto && c : arange(connectivity.getNbComponent())) {
          std::size_t node_tag;
          sstr_elem >> node_tag;

          AKANTU_DEBUG_ASSERT(node_tag >= file.first_node_number,
                              "Node number not in range : line "
                                  << file.current_line);

          AKANTU_DEBUG_ASSERT(node_tag <= file.last_node_number,
                              "Node number not in range : line "
                                  << file.current_line);

          node_tag = file.node_tags[node_tag];
          local_connect(read_order[c]) = node_tag;
        }
        connectivity.push_back(local_connect);
        elem.element = connectivity.size() - 1;
        file.element_tags[elem_tag] = elem;

        bool first = true;
        for (auto it = range.first; it != range.second; ++it) {
          auto phys_it = this->physical_names.find(it->second);
          if (first) {
            data0(elem.element) =
                it->second; // for compatibility with version 2
            if (phys_it != this->physical_names.end()) {
              physical_data(elem.element) = phys_it->second;
            }
            first = false;
          }
          if (phys_it != this->physical_names.end()) {
            file.mesh.getElementGroup(phys_it->second).add(elem, true, false);
          }
        }
      }
    }

    for (auto && element_group : file.mesh.iterateElementGroups()) {
      element_group.getNodeGroup().optimize();
    }
  };
}

/* -------------------------------------------------------------------------- */
void MeshIOMSH::read(const std::string & filename, Mesh & mesh) {

  File file(filename, mesh);

  std::map<std::string, std::function<void(const std::string &)>> readers;

  readers["$MeshFormat"] = [&](const std::string & /*unused*/) {
    auto && sstr = file.get_line();

    double version;
    int format;
    sstr >> version >> format;

    int major = std::trunc(version);
    int minor = std::round(10 * (version - major));
    file.version = major * 1000 + minor;

    if (format != 0) {
      AKANTU_ERROR("This reader can only read ASCII files.");
    }

    if (file.version > 2000) {
      sstr >> file.size_of_size_t;
      if (file.size_of_size_t > int(sizeof(UInt))) {
        AKANTU_DEBUG_INFO("The size of the indexes in akantu might be to small "
                          "to read this file (akantu "
                          << sizeof(UInt) << " vs. msh file "
                          << file.size_of_size_t << ")");
      }
    }

    if (file.version < 4000) {
      this->populateReaders2(file, readers);
    } else {
      this->populateReaders4(file, readers);
    }
  };

  auto && read_data = [&](auto && entity_tags, auto && get_data,
                          auto && read_data) {
    auto read_data_tags = [&](auto x) {
      UInt number_of_tags{0};
      file.read_line(number_of_tags);
      std::vector<decltype(x)> tags(number_of_tags);

      for (auto && tag : tags) {
        file.read_line(tag);
      }
      return tags;
    };

    auto && string_tags = read_data_tags(std::string{});
    auto && real_tags [[gnu::unused]] = read_data_tags(double{});
    auto && int_tags = read_data_tags(int{});

    for (auto & s : string_tags) {
      s = trim(s, '"');
    }

    auto id = string_tags[0];
    auto size = int_tags[2];
    auto nb_component = int_tags[1];
    auto & data = get_data(id, size, nb_component);

    for (auto n [[gnu::unused]] : arange(size)) {
      auto && sstr = file.get_line();

      size_t tag;
      sstr >> tag;
      const auto & entity = entity_tags[tag];
      read_data(entity, sstr, data, nb_component);
    }
  };

  readers["$NodeData"] = [&](const std::string & /*unused*/) {
    /* $NodeData
       numStringTags(ASCII int)
       stringTag(string) ...
       numRealTags(ASCII int)
       realTag(ASCII double) ...
       numIntegerTags(ASCII int)
       integerTag(ASCII int) ...
       nodeTag(size_t) value(double) ...
       ...
       $EndNodeData */
    read_data(
        file.node_tags,
        [&](auto && id, auto && size [[gnu::unused]],
            auto && nb_component [[gnu::unused]]) -> Array<double> & {
          auto & data =
              file.mesh.template getNodalData<double>(id, nb_component);
          data.resize(size);
          return data;
        },
        [&](auto && node, auto && sstr, auto && data,
            auto && nb_component [[gnu::unused]]) {
          for (auto c : arange(nb_component)) {
            sstr >> data(node, c);
          }
        });
  };

  readers["$ElementData"] = [&](const std::string & /*unused*/) {
    /* $ElementData
       numStringTags(ASCII int)
       stringTag(string) ...
       numRealTags(ASCII int)
       realTag(ASCII double) ...
       numIntegerTags(ASCII int)
       integerTag(ASCII int) ...
       elementTag(size_t) value(double) ...
       ...
       $EndElementData
    */
    read_data(
        file.element_tags,
        [&](auto && id, auto && size [[gnu::unused]],
            auto && nb_component
            [[gnu::unused]]) -> ElementTypeMapArray<double> & {
          file.mesh.template getElementalData<double>(id);
          return file.mesh.template getElementalData<double>(id);
        },
        [&](auto && element, auto && sstr, auto && data, auto && nb_component) {
          if (not data.exists(element.type)) {
            data.alloc(mesh.getNbElement(element.type), nb_component,
                       element.type, element.ghost_type);
          }
          auto & data_array = data(element.type);
          for (auto c : arange(nb_component)) {
            sstr >> data_array(element.element, c);
          }
        });
  };

  readers["$ElementNodeData"] = [&](const std::string & /*unused*/) {
    /* $ElementNodeData
       numStringTags(ASCII int)
       stringTag(string) ...
       numRealTags(ASCII int)
       realTag(ASCII double) ...
       numIntegerTags(ASCII int)
       integerTag(ASCII int) ...
       elementTag(size_t) value(double) ...
       ...
       $EndElementNodeData
    */
    read_data(
        file.element_tags,
        [&](auto && id, auto && size [[gnu::unused]],
            auto && nb_component
            [[gnu::unused]]) -> ElementTypeMapArray<double> & {
          file.mesh.template getElementalData<double>(id);
          auto & data = file.mesh.template getElementalData<double>(id);
          data.isNodal(true);
          return data;
        },
        [&](auto && element, auto && sstr, auto && data, auto && nb_component) {
          int nb_nodes_per_element;
          sstr >> nb_nodes_per_element;
          if (not data.exists(element.type)) {
            data.alloc(mesh.getNbElement(element.type),
                       nb_component * nb_nodes_per_element, element.type,
                       element.ghost_type);
          }
          auto & data_array = data(element.type);
          for (auto c : arange(nb_component)) {
            sstr >> data_array(element.element, c);
          }
        });
  };

  readers["$PhysicalNames"] = [&](const std::string & /*unused*/) {
    file.has_physical_names = true;
    int num_of_phys_names;
    file.read_line(num_of_phys_names); /// the format line

    for (auto k [[gnu::unused]] : arange(num_of_phys_names)) {
      int phys_name_id;
      int phys_dim;
      std::string phys_name;
      file.read_line(phys_dim, phys_name_id, std::quoted(phys_name));

      this->physical_names[phys_name_id] = phys_name;
    }
  };

  readers["Unsupported"] = [&](const std::string & _block) {
    std::string block = _block.substr(1);
    AKANTU_DEBUG_WARNING("Unsupported block_kind " << block << " at line "
                                                   << file.current_line);
    auto && end_block = "$End" + block;
    while (file.line != end_block) {
      file.get_line();
    }
  };

  while (file.good()) {
    std::string block;
    file.read_line(block);

    auto && it = readers.find(block);

    if (it != readers.end()) {
      it->second(block);

      std::string end_block;
      file.read_line(end_block);
      block = block.substr(1);

      if (end_block != "$End" + block) {
        AKANTU_EXCEPTION("The reader failed to properly read the block "
                         << block << ". Expected a $End" << block << " at line "
                         << file.current_line);
      }
    } else if (not block.empty()) {
      readers["Unsupported"](block);
    }
  }

  if (file.version < 4000) {
    this->constructPhysicalNames("tag_0", mesh);
    if (file.has_physical_names) {
      mesh.createGroupsFromMeshData<std::string>("physical_names");
    }
  }

  MeshUtils::fillElementToSubElementsData(mesh);
}

/* -------------------------------------------------------------------------- */
void MeshIOMSH::write(const std::string & filename, const Mesh & mesh) {
  std::ofstream outfile;
  const Array<Real> & nodes = mesh.getNodes();

  outfile.open(filename.c_str());

  outfile << "$MeshFormat"
          << "\n";
  outfile << "2.2 0 8"
          << "\n";
  outfile << "$EndMeshFormat"
          << "\n";

  outfile << std::setprecision(std::numeric_limits<Real>::digits10);
  outfile << "$Nodes"
          << "\n";

  outfile << nodes.size() << "\n";

  outfile << std::uppercase;
  for (Int i = 0; i < nodes.size(); ++i) {
    auto offset = i * nodes.getNbComponent();
    outfile << i + 1;
    for (Int j = 0; j < nodes.getNbComponent(); ++j) {
      outfile << " " << nodes.data()[offset + j];
    }

    for (Int p = nodes.getNbComponent(); p < 3; ++p) {
      outfile << " " << 0.;
    }
    outfile << "\n";
    ;
  }
  outfile << std::nouppercase;
  outfile << "$EndNodes"
          << "\n";

  outfile << "$Elements"
          << "\n";
  Int nb_elements = 0;
  for (auto && type :
       mesh.elementTypes(_all_dimensions, _not_ghost, _ek_not_defined)) {
    const auto & connectivity = mesh.getConnectivity(type, _not_ghost);
    nb_elements += connectivity.size();
  }
  outfile << nb_elements << "\n";

  std::map<Element, size_t> element_to_msh_element;

  Idx element_idx = 1;
  auto element = ElementNull;
  for (auto && type :
       mesh.elementTypes(_all_dimensions, _not_ghost, _ek_not_defined)) {
    const auto & connectivity = mesh.getConnectivity(type, _not_ghost);
    element.type = type;

    const Int * tag[2] = {nullptr, nullptr};
    if (mesh.hasData<Int>("tag_0", type, _not_ghost)) {
      const auto & data_tag_0 = mesh.getData<Int>("tag_0", type, _not_ghost);
      tag[0] = data_tag_0.data();
    }

    if (mesh.hasData<Int>("tag_1", type, _not_ghost)) {
      const auto & data_tag_1 = mesh.getData<Int>("tag_1", type, _not_ghost);
      tag[1] = data_tag_1.data();
    }

    for (auto && data :
         enumerate(make_view(connectivity, connectivity.getNbComponent()))) {
      element.element = std::get<0>(data);
      const auto & conn = std::get<1>(data);
      element_to_msh_element.insert(std::make_pair(element, element_idx));

      outfile << element_idx << " " << _akantu_to_msh_element_types[type]
              << " 2";

      /// \todo write the real data in the file
      for (Int t = 0; t < 2; ++t) {
        if (tag[t] != nullptr) {
          outfile << " " << tag[t][element.element];
        } else {
          outfile << " 0";
        }
      }

      for (auto && c : conn) {
        outfile << " " << c + 1;
      }
      outfile << "\n";
      element_idx++;
    }
  }
  outfile << "$EndElements"
          << "\n";

  if (mesh.hasData(MeshDataType::_nodal)) {
    auto && tags = mesh.getTagNames();
    for (auto && tag : tags) {
      auto type = mesh.getTypeCode(tag, MeshDataType::_nodal);
      if (type != MeshDataTypeCode::_real) {
        AKANTU_DEBUG_WARNING(
            "The field "
            << tag << " is ignored by the MSH writer, msh files do not support "
            << type << " data");
        continue;
      }
      auto && data = mesh.getNodalData<double>(tag);
      outfile << "$NodeData"
              << "\n";
      outfile << "1"
              << "\n";
      outfile << "\"" << tag << "\"\n";
      outfile << "1\n0.0"
              << "\n";
      outfile << "3\n0"
              << "\n";
      outfile << data.getNbComponent() << "\n";
      outfile << data.size() << "\n";
      for (auto && d : enumerate(make_view(data, data.getNbComponent()))) {
        outfile << std::get<0>(d) + 1;
        for (auto && v : std::get<1>(d)) {
          outfile << " " << v;
        }
        outfile << "\n";
      }
      outfile << "$EndNodeData"
              << "\n";
    }
  }

  if (mesh.hasData(MeshDataType::_elemental)) {
    auto && tags = mesh.getTagNames();
    for (auto && tag : tags) {
      auto && data = mesh.getElementalData<double>(tag);
      auto type = mesh.getTypeCode(tag, MeshDataType::_elemental);
      if (type != MeshDataTypeCode::_real) {
        AKANTU_DEBUG_WARNING(
            "The field "
            << tag << " is ignored by the MSH writer, msh files do not support "
            << type << " data");
        continue;
      }
      if (data.isNodal()) {
        continue;
      }

      auto size = data.size();
      if (size == 0) {
        continue;
      }
      auto && nb_components = data.getNbComponents();
      auto nb_component = nb_components(*(data.elementTypes().begin()));

      outfile << "$ElementData"
              << "\n";
      outfile << "1"
              << "\n";
      outfile << "\"" << tag << "\"\n";
      outfile << "1\n0.0"
              << "\n";
      outfile << "3\n0"
              << "\n";
      outfile << nb_component << "\n";
      outfile << size << "\n";

      Element element;
      for (auto type : data.elementTypes()) {
        element.type = type;
        for (auto && _ :
             enumerate(make_view(data(type), nb_components(type)))) {
          element.element = std::get<0>(_);
          outfile << element_to_msh_element[element];
          for (auto && v : std::get<1>(_)) {
            outfile << " " << v;
          }
          outfile << "\n";
        }
      }
      outfile << "$EndElementData"
              << "\n";
    }
  }

  outfile.close();
}

/* --------------------------------------------------------------------------
 */

} // namespace akantu
