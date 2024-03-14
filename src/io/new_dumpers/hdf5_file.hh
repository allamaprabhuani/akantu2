/**
 * @file   dumper_hdf5.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Oct 03 2017
 *
 * @brief Dump data in xdmf format
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
#include "aka_iterators.hh"
#include "dumper_field.hh"
#include "dumper_file_base.hh"
#include "support.hh"
/* -------------------------------------------------------------------------- */
#include <array>
#include <filesystem>
#include <hdf5.h>
#include <sstream>
#include <typeindex>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define ENDIANNESS(x) (x##LE)
#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
#define ENDIANNESS(x) (x##BE)
#else
#error "Does not know from which end to open the eggs"
#endif

namespace akantu {
namespace dumper {
  namespace {
    template <class T> constexpr std::type_index TIDX(T x) {
      return std::type_index(typeid(x));
    }

    hid_t datatype_id_out(const std::type_index & type_idx) {
      const static std::unordered_map<std::type_index, hid_t> type_ids{
          {TIDX(int{}), sizeof(int) == 4 ? ENDIANNESS(H5T_STD_I32)
                                         : ENDIANNESS(H5T_STD_I64)},
          {TIDX(unsigned{}), sizeof(int) == 4 ? ENDIANNESS(H5T_STD_U32)
                                              : ENDIANNESS(H5T_STD_U64)},
          {TIDX(int32_t{}), ENDIANNESS(H5T_STD_I32)},
          {TIDX(int64_t{}), ENDIANNESS(H5T_STD_I64)},
          {TIDX(uint32_t{}), ENDIANNESS(H5T_STD_U32)},
          {TIDX(uint64_t{}), ENDIANNESS(H5T_STD_U64)},
          {TIDX(bool{}), ENDIANNESS(H5T_STD_U8)},
          {TIDX(float{}), ENDIANNESS(H5T_IEEE_F32)},
          {TIDX(double{}), ENDIANNESS(H5T_IEEE_F64)}};

      return type_ids.at(type_idx);
    }

    hid_t datatype_id_in(const std::type_index & type_idx) {
      const static std::unordered_map<std::type_index, hid_t> type_ids{
          {TIDX(int{}), H5T_NATIVE_INT},
          {TIDX(unsigned{}), H5T_NATIVE_UINT},
          {TIDX(int32_t{}), H5T_NATIVE_INT32},
          {TIDX(int64_t{}), H5T_NATIVE_INT64},
          {TIDX(uint32_t{}), H5T_NATIVE_UINT32},
          {TIDX(uint64_t{}), H5T_NATIVE_UINT64},
          {TIDX(bool{}), H5T_NATIVE_CHAR},
          {TIDX(float{}), H5T_NATIVE_FLOAT},
          {TIDX(double{}), H5T_NATIVE_DOUBLE},
      };

      return type_ids.at(type_idx);
    }
  } // namespace

  namespace HDF5 {
    namespace fs = std::filesystem;

    enum class EntityType { _group, _dataset, _file, _dataspace, _link };

    struct Entity {
      hid_t id{};
      fs::path path;
      EntityType type;
      bool is_open{false};

      Entity(const fs::path & path, EntityType type) : path(path), type(type) {}

      Entity(const Entity & other) = delete;
      Entity(Entity && other) = delete;
      Entity & operator=(const Entity & other) = delete;
      Entity & operator=(Entity && other) = delete;

      void close() {
        if (not is_open) {
          return;
        }
        static const std::unordered_map<EntityType, std::function<void(hid_t)>>
            close_func{{EntityType::_group, H5Gclose},
                       {EntityType::_file, H5Fclose},
                       {EntityType::_dataset, H5Dclose},
                       {EntityType::_dataspace, H5Sclose}};
        if (auto it = close_func.find(type); it != close_func.end()) {
          it->second(id);
        }

        is_open = false;
      }

      void open(const fs::path & path, unsigned mode = H5F_ACC_RDWR,
                hid_t apl_id = H5P_DEFAULT) {
        if (is_open) {
          return;
        }
        switch (type) {
        case EntityType::_file:
          this->id = H5Fopen(path.c_str(), mode, apl_id);
          break;
        default:
          AKANTU_TO_IMPLEMENT();
        }
        is_open = true;
      }

      void open(hid_t group_id, hid_t apl_id = H5P_DEFAULT) {
        if (is_open) {
          return;
        }
        switch (type) {
        case EntityType::_group:
          this->id = H5Gopen(group_id, path.c_str(), apl_id);
          break;
        case EntityType::_dataset:
          this->id = H5Dopen(group_id, path.c_str(), apl_id);
          break;
        default:
          AKANTU_TO_IMPLEMENT();
        }
        is_open = true;
      }

      ~Entity() { close(); }
    };

    class File : public FileBase {
    public:
      File(SupportBase & support, const fs::path & path,
           hid_t fapl_id = H5P_DEFAULT)
          : FileBase(support) {
        auto path_wof = path;
        path_wof.remove_filename();
        if (not fs::exists(path_wof)) {
          fs::create_directories(path_wof);
        }

        auto file = std::make_unique<Entity>("/", EntityType::_file);
        file->id = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
        file->is_open = true;

        support.addProperty("hdf5_file", path.string());

        entities.push_back(std::move(file));

        /* Turn off error handling */
        H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
      }

      void close() { entities[0]->close(); }

      void open(const fs::path & path, hid_t fapl_id = H5P_DEFAULT) {
        entities[0]->open(path, H5F_ACC_RDWR, fapl_id);
      }

    protected:
      auto & openGroup(const std::string & path) {
        auto & group = *entities.back();

        auto && new_group =
            std::make_unique<Entity>(group.path, EntityType::_group);
        new_group->path /= path;

        auto status =
            H5Oexists_by_name(group.id, new_group->path.c_str(), H5P_DEFAULT);
        if (status <= 0) {
          AKANTU_DEBUG_INFO("DumperHDF5: creating group " << path);
          new_group->id = H5Gcreate(group.id, new_group->path.c_str(),
                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          new_group->is_open = true;
        } else {
          AKANTU_DEBUG_INFO("DumperHDF5: opening group " << path << " in "
                                                         << group.path);
          new_group->open(group.id);
        }
        entities.push_back(std::move(new_group));
        return *entities.back();
      }

      auto & createSymlink(const std::string & path, const std::string & link) {
        auto & group = *entities.back();

        auto && new_link =
            std::make_unique<Entity>(group.path, EntityType::_link);
        new_link->path /= path;

        auto status [[maybe_unused]] =
            H5Lcreate_soft(link.c_str(), group.id, new_link->path.c_str(),
                           H5P_DEFAULT, H5P_DEFAULT);

        AKANTU_DEBUG_ASSERT(status >= 0, "Could not create HDF5 link "
                                             << new_link->path.c_str() << " -> "
                                             << link.c_str()
                                             << ", status=" << status);

        entities.push_back(std::move(new_link));
        return *entities.back();
      }

    private:
      template <class T> bool unchanged(const T & t) {
        return (t.hasProperty("hdf5_release") and
                t.getRelease() ==
                    t.template getProperty<Release>("hdf5_release"));
      }

      void createElementalDataspace(FieldBase & field) {
        auto && field_element = aka::as_type<FieldElementalArrayBase>(field);

        auto && global_sizes =
            field.getSupport().getProperty<ElementTypeMap<hsize_t>>(
                "hdf5_element_global_sizes");
        std::array<hsize_t, 2> dims{
            hsize_t(global_sizes(field_element.getElementType())),
            field.getProperty<hsize_t>("nb_components")};

        std::array<hsize_t, 2> data_dims{field.getProperty<hsize_t>("size"),
                                         dims[1]};

        auto && offsets =
            field.getSupport().getProperty<ElementTypeMap<hsize_t>>(
                "hdf5_element_offsets");
        std::array<hsize_t, 2> offset{offsets(field_element.getElementType()),
                                      0};

        auto memoryspace_id =
            H5Screate_simple(data_dims.size(), data_dims.data(), nullptr);

        auto dataspace_id = H5Screate_simple(dims.size(), dims.data(), nullptr);

        auto filespace_id = H5Scopy(dataspace_id);
        H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset.data(),
                            nullptr, data_dims.data(), nullptr);

        field.addProperty<hid_t>("hdf5_memoryspace", memoryspace_id);
        field.addProperty<hid_t>("hdf5_filespace", filespace_id);
        field.addProperty<hid_t>("hdf5_dataspace", dataspace_id);
      }

      void createNodalDataspace(FieldBase & field) {
        std::array<hsize_t, 2> dims{
            field.getProperty<hsize_t>("size"),
            field.getProperty<hsize_t>("nb_components")};
        std::array<hsize_t, 2> data_dims{
            field.getSupport().getProperty<hsize_t>("hdf5_node_global_sizes"),
            dims[1]};

        auto memoryspace_id =
            H5Screate_simple(dims.size(), dims.data(), nullptr);
        auto dataspace_id =
            H5Screate_simple(data_dims.size(), data_dims.data(), nullptr);
        auto filespace_id = H5Scopy(dataspace_id);

        H5Sselect_none(memoryspace_id);
        H5Sselect_none(filespace_id);

        auto && slabs =
            field.getSupport().getProperty<PropertiesManager::slabs_type>(
                "hdf5_nodes_slabs");
        std::array<hsize_t, 2> slab_dims{0, dims[1]};
        std::array<hsize_t, 2> offsets{0, 0};

        for (auto && [local_offset, global_offset, length] : slabs) {
          offsets[0] = local_offset;
          slab_dims[0] = length;
          H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_OR, offsets.data(),
                              nullptr, slab_dims.data(), nullptr);

          offsets[0] = global_offset;
          H5Sselect_hyperslab(filespace_id, H5S_SELECT_OR, offsets.data(),
                              nullptr, slab_dims.data(), nullptr);
        }

        field.addProperty<hid_t>("hdf5_memoryspace", memoryspace_id);
        field.addProperty<hid_t>("hdf5_filespace", filespace_id);
        field.addProperty<hid_t>("hdf5_dataspace", dataspace_id);
      }

      void createDataspace(FieldBase & field) {
        field.addProperty<hid_t>("hdf5_filespace", H5I_INVALID_HID);

        using dumper::FieldType;
        switch (field.getFieldType()) {
        case FieldType::_node_array:
          createNodalDataspace(field);
          break;
        case FieldType::_element_array:
          createElementalDataspace(field);
          break;
        case FieldType::_element_map_array: /* FALLTHRU */
        case FieldType::_internal_field:    /* FALLTHRU */
        case FieldType::_not_defined:       /* FALLTHRU */
        default:
          AKANTU_EXCEPTION("The field type is not properly defined");
          break;
        }
      }

      auto createDataSet(FieldBase & field) {
        auto & group = *entities.back();
        auto data_set =
            std::make_unique<Entity>(group.path, EntityType::_dataset);
        data_set->path /= field.getName();

        if (not field.hasProperty("hdf5_dataspace") or
            not field.hasProperty("hdf5_description_release") or
            field.getProperty<Release>("hdf5_description_release") !=
                field.getSupport().getProperty<Release>(
                    "hdf5_description_release")) {
          createDataspace(field);

          field.addProperty("hdf5_description_release",
                            field.getSupport().getProperty<Release>(
                                "hdf5_description_release"));
        }

        auto status =
            H5Oexists_by_name(group.id, data_set->path.c_str(), H5P_DEFAULT);

        if (status <= 0) {

          AKANTU_DEBUG_INFO("DumperHDF5: creating data-set "
                            << field.getName() << " in "
                            << group.path.generic_string());
          data_set->id = H5Dcreate(group.id, data_set->path.c_str(),
                                   datatype_id_out(field.type()),
                                   field.getProperty<hid_t>("hdf5_dataspace"),
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          data_set->is_open = true;
        } else {
          AKANTU_DEBUG_INFO("HDF5: opening existing data-set "
                            << data_set->path);
          data_set->open(group.id);
        }

        field.addProperty("hdf5_path", data_set->path.generic_string());
        return data_set;
      }

      void createDataDescriptions(Support<Mesh> & support) {
        // Information needed for elements
        std::vector<hsize_t> offsets_v, global_sizes_v;
        for (auto type : support.elementTypes()) {
          global_sizes_v.push_back(support.getNbElements(type));
        }
        global_sizes_v.push_back(support.getNbLocalNodes());
        support.getCommunicator().exclusiveScan(global_sizes_v, offsets_v);
        support.getCommunicator().allReduce(global_sizes_v);

        ElementTypeMap<hsize_t> offsets, global_sizes;
        for (auto && [i, type] : enumerate(support.elementTypes())) {
          offsets(type) = offsets_v[i];
          global_sizes(type) = global_sizes_v[i];
        }

        const auto & mesh = support.getMesh();

        // Information needed for nodes
        PropertiesManager::slabs_type slabs;
        hsize_t nb_nodes = support.getNbNodes();
        hsize_t node = 0;
        while (node < nb_nodes) {
          auto first_node = node;
          auto global_counter = mesh.getNodeGlobalId(node);
          hsize_t first_global_node = global_counter;
          while (node < nb_nodes and
                 global_counter == mesh.getNodeGlobalId(node) and
                 mesh.isLocalOrMasterNode(node)) {
            ++node;
            ++global_counter;
          }

          if (first_node != node) {
            slabs.push_back({first_node, first_global_node, node - first_node});
          }

          while (node < nb_nodes and not mesh.isLocalOrMasterNode(node)) {
            ++node;
          }
        }

        support.addProperty("hdf5_element_offsets", offsets);
        support.addProperty("hdf5_element_global_sizes", global_sizes);
        support.addProperty("hdf5_node_global_sizes", global_sizes_v.back());
        support.addProperty("hdf5_nodes_slabs", slabs);
        support.addProperty("hdf5_description_release", support.getRelease());
      }

    protected:
      void dump(FieldArrayBase & field) override {
        field.addProperty("size", field.size());
        field.addProperty("nb_components", field.getNbComponent());

        auto dataset = createDataSet(field);

#if defined(AKANTU_USE_MPI)
        auto dxpl_id = H5Pcreate(H5P_DATASET_XFER);
        //        H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
        H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_INDEPENDENT);
#else
        auto dxpl_id = H5P_DEFAULT
#endif

        AKANTU_DEBUG_INFO("HDF5: writing dataset " << dataset->path);
        H5Dwrite(dataset->id, datatype_id_in(field.type()),
                 field.getProperty<hid_t>("hdf5_memoryspace"),
                 field.getProperty<hid_t>("hdf5_filespace"), dxpl_id,
                 field.data());

#if defined(AKANTU_USE_MPI)
        H5Pclose(dxpl_id);
#endif

        field.addProperty("hdf5_release", field.getRelease());
      }

      void dump(FieldElementMapArrayBase & field) override {
        if (unchanged(field)) {
          createSymlink(field.getName(),
                        field.getProperty<std::string>("hdf5_path"));
          entities.pop_back();
          return;
        }

        auto && group = openGroup(field.getName());
        field.addProperty("hdf5_path", group.path.generic_string());

        for (auto && type : field.elementTypes(_not_ghost)) {
          auto & array = field.array(type);
          if (not array.hasProperty("name")) {
            array.addProperty("name",
                              field.getName() + ":" + std::to_string(type));
          }
          dump(array);
        }
        field.addProperty("hdf5_release", field.getRelease());
        entities.pop_back();
      }

      void dump(Support<Mesh> & support) override {
        openGroup(support.getName());

        if (unchanged(support)) {
          createSymlink("topology",
                        support.getProperty<std::string>("hdf5_path"));
          entities.pop_back();
        } else {
          createDataDescriptions(support);

          auto & topology = openGroup("topology");
          support.addProperty("hdf5_path", topology.path.generic_string());

          openGroup("nodes");
          auto && nodes = support.getNodes();
          dump(nodes);
          entities.pop_back();

          auto && connectivities = support.getConnectivities();
          dump(connectivities);
          entities.pop_back();

          support.addProperty("hdf5_release", support.getRelease());
        }

        // auto && data_grp = Group(loc_id, path_ + "/data");

        // openGroup("groups");

        // for (auto && [_, support] : support.getSubSupports()) {
        //   dump(*support);
        // }

        // // close groups
        // entities.pop_back();

        // for (auto && [_, field] : support.getFields()) {
        //   FileBase::dump(*field);
        // }

        // close mesh
        entities.pop_back();
      }

      void dump(Support<ElementGroup> & support) override {
        const auto & parent_support = support.getParentSupport();
        openGroup(support.getName());

        if (unchanged(support)) {
          createSymlink("topology",
                        support.getProperty<std::string>("hdf5_path"));
          entities.pop_back();
        } else {
          support.addProperty(
              "hdf5_file",
              parent_support.getProperty<std::string>("hdf5_file"));

          auto & topology = openGroup("topology");
          support.addProperty("hdf5_path", topology.path.generic_string());
          dump(support.getElements());
          dump(support.getConnectivities());
          entities.pop_back();
          support.addProperty("hdf5_release", support.getRelease());
        }

        for (auto && [_, field] : support.getFields()) {
          FileBase::dump(*field);
        }

        entities.pop_back();
      }

      void dump(SupportBase & support) override { FileBase::dump(support); }

      void dump() override {
        openGroup("steps");

        if (support.hasProperty("time")) {
          openGroup(std::to_string(support.getProperty<double>("time")));
        }

        FileBase::dump(support);

        if (support.hasProperty("time")) {
          entities.pop_back();
        }
        entities.pop_back();
      }

    private:
      std::vector<std::unique_ptr<Entity>> entities;
    };

  } // namespace HDF5
} // namespace dumper
} // namespace akantu
