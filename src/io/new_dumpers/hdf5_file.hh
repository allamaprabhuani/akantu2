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
#define TIDX(x) std::type_index(typeid(x))
    hid_t datatype_id_out(const std::type_index & type_idx) {
      const static std::unordered_map<std::type_index, hid_t> type_ids{
          {TIDX(int), sizeof(int) == 4 ? ENDIANNESS(H5T_STD_I32)
                                       : ENDIANNESS(H5T_STD_I64)},
          {TIDX(unsigned int), sizeof(int) == 4 ? ENDIANNESS(H5T_STD_U32)
                                                : ENDIANNESS(H5T_STD_U64)},
          {TIDX(int32_t), ENDIANNESS(H5T_STD_I32)},
          {TIDX(int64_t), ENDIANNESS(H5T_STD_I64)},
          {TIDX(uint32_t), ENDIANNESS(H5T_STD_U32)},
          {TIDX(uint64_t), ENDIANNESS(H5T_STD_U64)},
          {TIDX(bool), ENDIANNESS(H5T_STD_U8)},
          {TIDX(float), ENDIANNESS(H5T_IEEE_F32)},
          {TIDX(double), ENDIANNESS(H5T_IEEE_F64)},
      };

      return type_ids.at(type_idx);
    }

    hid_t datatype_id_in(const std::type_index & type_idx) {
      const static std::unordered_map<std::type_index, hid_t> type_ids{
          {TIDX(int), H5T_NATIVE_INT},
          {TIDX(unsigned int), H5T_NATIVE_UINT},
          {TIDX(int32_t), H5T_NATIVE_INT32},
          {TIDX(int64_t), H5T_NATIVE_INT64},
          {TIDX(uint32_t), H5T_NATIVE_UINT32},
          {TIDX(uint64_t), H5T_NATIVE_UINT64},
          {TIDX(bool), H5T_NATIVE_CHAR},
          {TIDX(float), H5T_NATIVE_FLOAT},
          {TIDX(double), H5T_NATIVE_DOUBLE},
      };

      return type_ids.at(type_idx);
    }
#undef TIDX
  } // namespace

  namespace HDF5 {
    namespace fs = std::filesystem;

    enum class EntityType { _group, _dataset, _file, _dataspace, _link };

    struct Entity {
      hid_t id{};
      fs::path path;
      EntityType type;

      Entity(const fs::path & path, EntityType type) : path(path), type(type) {}

      Entity(const Entity & other) = delete;
      Entity(Entity && other) = delete;
      Entity & operator=(const Entity & other) = delete;
      Entity & operator=(Entity && other) = delete;

      ~Entity() {
        switch (type) {
        case EntityType::_group:
          H5Gclose(id);
          break;
        case EntityType::_file:
          H5Fclose(id);
          break;
        case EntityType::_dataset:
          H5Dclose(id);
          break;
        case EntityType::_dataspace:
          H5Sclose(id);
          break;
        case EntityType::_link:
          break;
        }
      }
    };

    class File : public FileBase {
    public:
      File(SupportBase & support, const fs::path & path) : FileBase(support) {
        auto path_wof = path;
        path_wof.remove_filename();
        if (not fs::exists(path_wof)) {
          fs::create_directories(path_wof);
        }

        auto file = std::make_unique<Entity>("/", EntityType::_file);
        file->id =
            H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        support.addProperty("hdf5_file", path.string());

        entities.push_back(std::move(file));

        /* Turn off error handling */
        H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
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
        } else {
          AKANTU_DEBUG_INFO("DumperHDF5: opening group " << path << " in "
                                                         << group.path);
          new_group->id =
              H5Gopen(group.id, new_group->path.c_str(), H5P_DEFAULT);
        }
        entities.push_back(std::move(new_group));
        return *entities.back();
      }

      auto & createSymlink(const std::string & path, const std::string & link) {
        auto & group = *entities.back();

        auto && new_link =
            std::make_unique<Entity>(group.path, EntityType::_link);
        new_link->path /= path;

        auto status =
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
      auto createDataSet(FieldBase & field) {
        auto & group = *entities.back();
        auto data_set =
            std::make_unique<Entity>(group.path, EntityType::_dataset);
        data_set->path /= field.getName();
        auto status =
            H5Oexists_by_name(group.id, data_set->path.c_str(), H5P_DEFAULT);

        decltype(data_set) data_space;

        if (status <= 0) {
          std::array<hsize_t, 2> dims{
              field.getProperty<hsize_t>("size"),
              field.getProperty<hsize_t>("nb_components")};
          data_space = std::make_unique<Entity>("", EntityType::_dataspace);
          data_space->id = H5Screate_simple(dims.size(), dims.data(), nullptr);

          AKANTU_DEBUG_INFO("DumperHDF5: creating data-set "
                            << field.getName() << " in "
                            << group.path.generic_string());
          data_set->id = H5Dcreate(
              group.id, data_set->path.c_str(), datatype_id_out(field.type()),
              data_space->id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        } else {
          AKANTU_DEBUG_INFO("HDF5: opening existing data-set "
                            << data_set->path);
          data_set->id = H5Dopen(group.id, data_set->path.c_str(), H5P_DEFAULT);
        }

        field.addProperty("hdf5_path", data_set->path.generic_string());
        return data_set;
      }

    protected:
      void dump(FieldNodeArrayBase & field) override {
        field.addProperty("size", field.size());
        field.addProperty("nb_components", field.getNbComponent());

        auto data_set = createDataSet(field);

        AKANTU_DEBUG_INFO("HDF5: writing dataset " << data_set->path);
        H5Dwrite(data_set->id, datatype_id_in(field.type()), H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, field.data());
      }

      void dump(FieldElementMapArrayBase & field) override {
        Int nb_components = -1;
        auto sizes = field.size();
        auto size = sizes.first;
        for (auto && type : field.elementTypes()) {
          auto nb_new_components = field.getNbComponent(type);
          if (nb_new_components != nb_components and nb_components != -1) {
            nb_components = 1;
            size = sizes.second;
            break;
          }
          nb_components = nb_new_components;
        }

        field.addProperty("size", size);
        field.addProperty("nb_components", nb_components);

        auto data_set = createDataSet(field);

        for (auto && type : field.elementTypes()) {
          auto & array = field.array(type);

          // TODO: create correspoondig dataspaces

          AKANTU_DEBUG_INFO("HDF5: writing dataset " << data_set->path);
          H5Dwrite(data_set->id, datatype_id_in(field.type()), H5S_ALL, H5S_ALL,
                   H5P_DEFAULT, array.data());
        }
      }

      // void dump(FieldConnectivity & field) {
      //   auto && [total_element, total_size] = field.size();
      //   total_size += total_element;

      //   field.addProperty("nb_total_element", total_element);

      //   Array<UInt> total_array(total_size);

      //   UInt i{0};
      //   for (auto && type : field.elementTypes()) {
      //     auto && array = field.arrayTyped(type).getArray();
      //     for (auto && connectivity :
      //          make_view(array, array.getNbComponent())) {
      //       total_array[i++] = aka_type_to_dumper_type.at(type);
      //       for (auto && c : connectivity) {
      //         total_array[i++] = c;
      //       }
      //     }
      //   }

      //   auto && connectivity_field = make_field(total_array,
      //   field.getSupport()); connectivity_field->addProperty("name",
      //                                   field.getProperty<std::string>("name"));
      //   dump(*connectivity_field);

      //   field.addProperty<int>("size", total_size);
      //   field.addProperty<int>("nb_components", 1);

      //   field.addProperty(
      //       "hdf5_path",
      //       connectivity_field->getProperty<std::string>("hdf5_path"));

      //   // entities.pop_back();
      // }

      void dump(Support<Mesh> & support) override {
        if (support.hasProperty("hdf5_release") and
            support.getRelease() == support.getProperty<Int>("hdf5_release")) {
          createSymlink(support.getName(),
                        support.getProperty<std::string>("hdf5_path"));
          entities.pop_back();
        }

        auto & group = openGroup(support.getName());
        support.addProperty("hdf5_path", group.path.generic_string());

        openGroup("nodes");
        auto && nodes = support.getNodes();
        dump(nodes);
        entities.pop_back();

        openGroup("elements");
        auto && connectivities = support.getConnectivities();
        dump(connectivities);
        entities.pop_back();

        support.addProperty("hdf5_release", support.getRelease());

        // auto && data_grp = Group(loc_id, path_ + "/data");

        openGroup("groups");

        for (auto && support : support.getSubSupports()) {
          dump(*support.second);
        }

        // close groups
        entities.pop_back();

        for (auto && field : support.getFields()) {
          FileBase::dump(*field.second);
        }

        // close mesh
        entities.pop_back();
      }

      void dump(Support<ElementGroup> & support) override {
        const auto & parent_support = support.getParentSupport();
        auto & group = openGroup(support.getName());
        support.addProperty("hdf5_path", group.path.generic_string());
        support.addProperty(
            "hdf5_file", parent_support.getProperty<std::string>("hdf5_file"));
        dump(support.getElements());

        for (auto && field : support.getFields()) {
          FileBase::dump(*field.second);
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
