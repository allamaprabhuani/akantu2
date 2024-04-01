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

    /* ---------------------------------------------------------------------- */
    // TODO: make a compute field that transfer an arry in to and array out with
    // FieldFunction being a special case
    template <class T>
    class HDF5NodeFieldTemplateIterator
        : public FieldComputeArray<T, T, FieldNodalArrayBase> {
      using parent = FieldComputeArray<T, T, FieldNodalArrayBase>;

    public:
      HDF5NodeFieldTemplateIterator(
          const std::shared_ptr<
              FieldArrayTemplateBase<T, FieldNodalArrayBase>> & array_in,
          const SupportBase & support)
          : parent(array_in, support) {
        for (auto prop : this->array_in->getPropertiesList()) {
          this->addProperty(prop, this->array_in->getPropertyVariant(prop));
        }
      }

      ~HDF5NodeFieldTemplateIterator() {
        for (auto prop : this->getPropertiesList()) {
          this->array_in->addProperty(prop, this->getPropertyVariant(prop));
        }
      }

      [[nodiscard]] bool unchanged() const override {
        return parent::unchanged() and
               (not this->array_in->hasProperty("hdf5_description_release") or
                this->array_in->template getProperty<Release>(
                    "hdf5_description_release") == description_release);
      }

      [[nodiscard]] Int size() const override {
        return dynamic_cast<const SupportElements &>(
                   this->array_in->getSupport())
            .getNbLocalNodes();
      }

      void compute() override {
        this->array_in->update();
        auto && slabs =
            this->array_in->getSupport()
                .template getProperty<PropertiesManager::slabs_type>(
                    "hdf5_nodes_slabs");
        auto nb_component = this->getNbComponent();

        auto out_it = make_view(*this->array_out, nb_component).begin();
        auto begin_in_it =
            make_view(this->array_in->getArray(), nb_component).begin();
        for (auto && [lid, gid, length] : slabs) {
          auto in_it = begin_in_it + lid;
          for (auto _ [[maybe_unused]] : arange(length)) {
            *out_it = *in_it;
            ++out_it;
            ++in_it;
          }
        }

        if (this->array_in->hasProperty("hdf5_description_release")) {
          description_release = this->array_in->template getProperty<Release>(
              "hdf5_description_release");
        }
      }

    protected:
      Release description_release;
    };

    auto make_hdf5_field(FieldNodalArrayBase & field)
        -> std::unique_ptr<FieldNodalArrayBase> {
      std::unique_ptr<FieldNodalArrayBase> out_field;
      if (field.type() == TIDX(Int{})) {
        out_field = std::make_unique<HDF5NodeFieldTemplateIterator<Int>>(
            dynamic_cast<FieldArrayTemplateBase<Int, FieldNodalArrayBase> &>(
                field)
                .getSharedPointer(),
            field.getSupport());
      } else if (field.type() == TIDX(Real{})) {
        out_field = std::make_unique<HDF5NodeFieldTemplateIterator<Real>>(
            dynamic_cast<FieldArrayTemplateBase<Real, FieldNodalArrayBase> &>(
                field)
                .getSharedPointer(),
            field.getSupport());
      } else {
        AKANTU_EXCEPTION("Not implemented for type: "
                         << debug::demangle(field.type().name()));
      }

      out_field->addProperty("name", field.getProperty<ID>("name"));
      return out_field;
    }

    /* ---------------------------------------------------------------------- */
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
        // H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
      }

      ~File() { close(); }

      void close() {
        for (auto && entity : entities) {
          entity->close();
        }
      }

      void open(const fs::path & path, hid_t fapl_id = H5P_DEFAULT) {
        entities[0]->open(path, H5F_ACC_RDWR, fapl_id);
      }

    protected:
      auto & openGroup(const std::string & path) {
        auto & group = *entities.back();

        auto && new_group =
            std::make_unique<Entity>(group.path, EntityType::_group);
        new_group->path /= path;

        auto status = H5Lexists(group.id, new_group->path.c_str(), H5P_DEFAULT);
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

      /* ------------------------------------------------------------------ */
      decltype(auto) getSlabs(Support<Mesh> & support) {
        PropertiesManager::slabs_type slabs;
        const auto & mesh = support.getMesh();
        // Information needed for nodes
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

        for (auto && [local_offset, global_offset, length] : slabs) {
          std::cout << Communicator::getWorldCommunicator().whoAmI()
                    << " - local: " << local_offset
                    << " - global: " << global_offset << " - length: " << length
                    << std::endl;
        }

        std::sort(slabs.begin(), slabs.end(),
                  [](const auto & a, const auto & b) {
                    return std::get<1>(a) < std::get<1>(b);
                  });

        return slabs;
      }

      /* ------------------------------------------------------------------ */
      decltype(auto) getSlabs(Support<ElementGroup> & support) {
        PropertiesManager::slabs_type slabs;
        const auto & mesh = support.getMesh();
        // Information needed for nodes
        hsize_t nb_nodes = support.getNbNodes();
        hsize_t node = 0;
        while (node < nb_nodes) {
          auto first_node = node;
          while (node < nb_nodes and mesh.isLocalOrMasterNode(node)) {
            ++node;
          }

          if (first_node != node) {
            auto offset = support.getNodesOffsets();
            slabs.push_back(
                {first_node, offset + first_node, node - first_node});
          }

          while (node < nb_nodes and not mesh.isLocalOrMasterNode(node)) {
            ++node;
          }
        }

        return slabs;
      }

      /* -------------------------------------------------------------------- */
      void createDataDescriptions(SupportElements & support) {
        auto & support_base = dynamic_cast<SupportBase &>(support);

        // Information needed for elements
        if (not support_base.hasProperty("hdf5_description_release") or
            support_base.getProperty<Release>("hdf5_description_release") !=
                support_base.getRelease()) {
          support.updateOffsets();
        }

        PropertiesManager::slabs_type slabs;
        switch (support_base.getType()) {
        case SupportType::_mesh:
          slabs = getSlabs(aka::as_type<Support<Mesh>>(support));
          break;
        case SupportType::_element_group:
          slabs = getSlabs(aka::as_type<Support<ElementGroup>>(support));
          break;
        }

        support_base.addProperty("hdf5_nodes_slabs", slabs);
        support_base.addProperty("hdf5_description_release",
                                 support_base.getRelease());
      }

      /* -------------------------------------------------------------------- */
      void createElementalDataspace(FieldBase & field) {
        auto && field_element = aka::as_type<FieldElementalArrayBase>(field);

        auto && support =
            dynamic_cast<const SupportElements &>(field.getSupport());
        auto element_type = field_element.getElementType();
        auto ghost_type = field_element.getGhostType();

        hsize_t nb_global_elements =
            support.getNbGlobalElements(element_type, ghost_type);
        std::array<hsize_t, 2> data_dims{
            field.getProperty<hsize_t>("size"),
            field.getProperty<hsize_t>("nb_components")};

        std::array<hsize_t, 2> dims{nb_global_elements, data_dims[1]};
        field.addProperty("global_size", nb_global_elements);

        hsize_t element_offset =
            support.getElementsOffsets(element_type, ghost_type);
        std::array<hsize_t, 2> offset{element_offset, 0};

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
        auto && support =
            dynamic_cast<const SupportElements &>(field.getSupport());

        hsize_t nb_global_nodes = support.getNbGlobalNodes();

        std::array<hsize_t, 2> data_dims{
            field.getProperty<hsize_t>("size"),
            field.getProperty<hsize_t>("nb_components")};
        std::array<hsize_t, 2> dims{nb_global_nodes, data_dims[1]};

        field.addProperty("global_size", nb_global_nodes);

        auto memoryspace_id =
            H5Screate_simple(data_dims.size(), data_dims.data(), nullptr);
        auto dataspace_id = H5Screate_simple(dims.size(), dims.data(), nullptr);
        auto filespace_id = H5Scopy(dataspace_id);

        auto && slabs =
            field.getSupport().getProperty<PropertiesManager::slabs_type>(
                "hdf5_nodes_slabs");
        std::array<hsize_t, 2> slab_dims{0, dims[1]};
        std::array<hsize_t, 2> offsets{0, 0};

        if (slabs.empty()) {
          // H5Sselect_none(memoryspace_id);
          H5Sselect_none(filespace_id);
        }

        H5S_seloper_t select_op = H5S_SELECT_SET;
        for (auto && [local_offset, global_offset, length] : slabs) {
          // offsets[0] = local_offset;
          slab_dims[0] = length;
          // H5Sselect_hyperslab(memoryspace_id, select_op, offsets.data(),
          //                     nullptr, slab_dims.data(), nullptr);

          // AKANTU_DEBUG_ASSERT(offsets[0] + slab_dims[0] <= data_dims[0],
          //                     "The selected hyperslab does not fit in
          //                     memory");

          offsets[0] = global_offset;
          AKANTU_DEBUG_ASSERT(offsets[0] + slab_dims[0] <= dims[0],
                              "The selected hyperslab does not fit in file");
          H5Sselect_hyperslab(filespace_id, select_op, offsets.data(), nullptr,
                              slab_dims.data(), nullptr);

          select_op = H5S_SELECT_OR;
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

        createDataspace(field);

        auto status = H5Lexists(group.id, data_set->path.c_str(), H5P_DEFAULT);

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

    protected:
      void dump(FieldNodalArrayBase & field) override {
        auto hdf5_field = make_hdf5_field(field);
        dump(dynamic_cast<FieldArrayBase &>(*hdf5_field));
      }

      void dump(FieldArrayBase & field) override {
        field.addProperty("size", field.size());
        field.addProperty("nb_components", field.getNbComponent());

        auto dataset = createDataSet(field);

#if defined(AKANTU_USE_MPI)
        auto dxpl_id = H5Pcreate(H5P_DATASET_XFER);
        //        H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
        H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_INDEPENDENT);
#else
        auto dxpl_id = H5P_DEFAULT;
#endif

        auto memspace_id = field.getProperty<hid_t>("hdf5_memoryspace");
        auto filespace_id = field.getProperty<hid_t>("hdf5_filespace");
        auto dataspace_id = field.getProperty<hid_t>("hdf5_dataspace");

        field.update();

        AKANTU_DEBUG_INFO("HDF5: writing dataset " << dataset->path);
        H5Dwrite(dataset->id, datatype_id_in(field.type()), memspace_id,
                 filespace_id, dxpl_id,
                 field.size() == 0 ? nullptr : field.data());

#if defined(AKANTU_USE_MPI)
        H5Pclose(dxpl_id);
#endif

        H5Sclose(memspace_id);
        H5Sclose(filespace_id);
        H5Sclose(dataspace_id);

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

        openGroup("groups");
        for (auto && [_, support] : support.getSubSupports()) {
          createDataDescriptions(aka::as_type<SupportElements>(*support));
          dump(*support);
        }

        // close groups
        entities.pop_back();

        openGroup("data");
        for (auto && [_, field] : support.getFields()) {
          FileBase::dump(*field);
        }
        entities.pop_back();

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

        openGroup("data");
        for (auto && [_, field] : support.getFields()) {
          FileBase::dump(*field);
        }
        entities.pop_back();

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
