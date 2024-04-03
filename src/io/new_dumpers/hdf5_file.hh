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
#include "hdf5_entities.hh"
#include "support.hh"
/* -------------------------------------------------------------------------- */
#include <array>
#include <filesystem>
#include <sstream>
#include <typeindex>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

namespace akantu {
namespace dumper {
  namespace HDF5 {
    namespace fs = std::filesystem;

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
    static herr_t hdf5_error_handler(hid_t estack_id, void * client_data);
  } // namespace HDF5

  namespace fs = std::filesystem;
  /* ------------------------------------------------------------------------ */
  class H5File : public FileBase {
  protected:
    fs::path filepath_;

  public:
    H5File(SupportBase & support, const fs::path & path,
           hid_t fapl_id = H5P_DEFAULT)
        : FileBase(support), filepath_(path) {
      auto path_wof = path;
      path_wof.remove_filename();
      if (not fs::exists(path_wof)) {
        fs::create_directories(path_wof);
      }

      H5Eset_auto(H5E_DEFAULT, HDF5::hdf5_error_handler, this);

      auto file = std::make_unique<HDF5::File>(path);
      file->create(H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);

      support.addProperty("hdf5_file", path.string());

      entities.push_back(std::move(file));
    }

    ~H5File() { close(); }

    void close() {
      AKANTU_DEBUG_ASSERT(entities.size() == 1,
                          "Not all object where closed, or too many");

      entities[0]->close();
    }

    void open(hid_t fapl_id = H5P_DEFAULT) {
      aka::as_type<HDF5::File>(*entities[0]).open(H5F_ACC_RDWR, fapl_id);
    }

    auto filepath() const { return filepath_; }

  protected:
    auto & openGroup(const std::string & path) {
      auto & group = *entities.back();

      auto && new_group = std::make_unique<HDF5::Group>(path, group);
      new_group->createOrOpen();

      entities.push_back(std::move(new_group));
      return *entities.back();
    }

    auto & createSymlink(const std::string & path, const std::string & link) {
      auto & group = *entities.back();

      auto && new_link =
          std::make_unique<HDF5::Entity<HDF5::EntityType::_link>>(group.path /
                                                                  path);

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

      std::sort(slabs.begin(), slabs.end(), [](const auto & a, const auto & b) {
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
          slabs.push_back({first_node, offset + first_node, node - first_node});
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
    decltype(auto) createElementalDataspace(FieldBase & field,
                                            HDF5::Dataset & dataset) {
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

      dataset.memoryspace = std::make_unique<HDF5::Dataspace>(data_dims);
      dataset.dataspace = std::make_unique<HDF5::Dataspace>(dims);
      dataset.filespace = std::make_unique<HDF5::Dataspace>(*dataset.dataspace);

      H5Sselect_hyperslab(dataset.filespace->id, H5S_SELECT_SET, offset.data(),
                          nullptr, data_dims.data(), nullptr);
    }

    void createNodalDataspace(FieldBase & field, HDF5::Dataset & dataset) {
      auto && support =
          dynamic_cast<const SupportElements &>(field.getSupport());

      hsize_t nb_global_nodes = support.getNbGlobalNodes();

      std::array<hsize_t, 2> data_dims{
          field.getProperty<hsize_t>("size"),
          field.getProperty<hsize_t>("nb_components")};
      std::array<hsize_t, 2> dims{nb_global_nodes, data_dims[1]};

      field.addProperty("global_size", nb_global_nodes);

      dataset.memoryspace = std::make_unique<HDF5::Dataspace>(data_dims);
      dataset.dataspace = std::make_unique<HDF5::Dataspace>(dims);
      dataset.filespace = std::make_unique<HDF5::Dataspace>(*dataset.dataspace);

      auto && slabs =
          field.getSupport().getProperty<PropertiesManager::slabs_type>(
              "hdf5_nodes_slabs");
      std::array<hsize_t, 2> slab_dims{0, dims[1]};
      std::array<hsize_t, 2> offsets{0, 0};

      if (slabs.empty()) {
        H5Sselect_none(dataset.filespace->id);
      }

      H5S_seloper_t select_op = H5S_SELECT_SET;
      for (auto && [local_offset, global_offset, length] : slabs) {
        slab_dims[0] = length;

        offsets[0] = global_offset;
        AKANTU_DEBUG_ASSERT(offsets[0] + slab_dims[0] <= dims[0],
                            "The selected hyperslab does not fit in file");
        H5Sselect_hyperslab(dataset.filespace->id, select_op, offsets.data(),
                            nullptr, slab_dims.data(), nullptr);

        select_op = H5S_SELECT_OR;
      }
    }

    void createDataspace(FieldBase & field, HDF5::Dataset & dataset) {
      field.addProperty<hid_t>("hdf5_filespace", H5I_INVALID_HID);

      using dumper::FieldType;
      switch (field.getFieldType()) {
      case FieldType::_node_array:
        createNodalDataspace(field, dataset);
        break;
      case FieldType::_element_array:
        createElementalDataspace(field, dataset);
        break;
      default:
        AKANTU_EXCEPTION("The field type is not properly defined");
        break;
      }
    }

    auto & createDataSet(FieldBase & field) {
      auto & group = aka::as_type<HDF5::Group>(*entities.back());
      auto data_set =
          std::make_unique<HDF5::Dataset>(field.getName(), group, field.type());

      createDataspace(field, *data_set);

      data_set->createOrOpen();

      field.addProperty("hdf5_path", data_set->path.generic_string());

      entities.push_back(std::move(data_set));
      return aka::as_type<HDF5::Dataset>(*entities.back());
    }

    auto getDataTransferPropertyList() {
      auto dxpl = std::make_unique<HDF5::PropertyList>();
#if defined(AKANTU_USE_MPI)
      dxpl->create(H5P_DATASET_XFER);
      //        H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
      dxpl->setDxplMPIIO(H5FD_MPIO_INDEPENDENT);
#endif
      return dxpl;
    }

  protected:
    void dump(FieldNodalArrayBase & field) override {
      auto hdf5_field = HDF5::make_hdf5_field(field);
      dump(dynamic_cast<FieldArrayBase &>(*hdf5_field));
    }

    void dump(FieldArrayBase & field) override {
      field.addProperty("size", field.size());
      field.addProperty("nb_components", field.getNbComponent());

      auto & dataset = createDataSet(field);
      auto dxpl = getDataTransferPropertyList();

      field.update();

      dataset.write(field.size() == 0 ? nullptr : field.data(), *dxpl);

      field.addProperty("hdf5_release", field.getRelease());

      entities.pop_back(); // dataset
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
            "hdf5_file", parent_support.getProperty<std::string>("hdf5_file"));

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
    std::vector<std::unique_ptr<HDF5::EntityBase>> entities;
  };

  namespace HDF5 {
    static herr_t hdf5_error_walk(unsigned int n, const H5E_error2_t * err_desc,
                                  void * client_data) {
      const int MSG_SIZE = 256;
      auto & messages =
          *reinterpret_cast<std::vector<std::string> *>(client_data);
      char maj[MSG_SIZE];
      char min[MSG_SIZE];
      char cls[MSG_SIZE];

      if (H5Eget_class_name(err_desc->cls_id, cls, MSG_SIZE) < 0) {
        return -1;
      }

      if (H5Eget_msg(err_desc->maj_num, NULL, maj, MSG_SIZE) < 0) {
        return -1;
      }

      if (H5Eget_msg(err_desc->min_num, NULL, min, MSG_SIZE) < 0) {
        return -1;
      }

      messages.push_back(
          std::string(cls) + " [error: " + std::to_string(int(n)) + "]" +
          " - " + std::string(maj) + " (" + std::string(min) + ")" +
          (err_desc->func_name ? (std::string(err_desc->func_name) + " [" +
                                  std::string(err_desc->file_name) + ":" +
                                  std::to_string(err_desc->line) + "]")
                               : ""));
      return 0;
    };

    static herr_t hdf5_error_handler(hid_t estack_id, void * client_data) {
      auto & file = *reinterpret_cast<H5File *>(client_data);

      std::vector<std::string> messages;
      H5Ewalk(estack_id, H5E_WALK_DOWNWARD, hdf5_error_walk, &messages);

      std::string errors;
      for (auto && [n, msg] : enumerate(messages)) {
        errors += "\n";
        errors += std::string(" ", AKANTU_INDENT) + msg;
      }

      AKANTU_EXCEPTION("HDF5 errors while accessing file "
                       << file.filepath().generic_string() << ": " << errors);

      return 0;
    }
  } // namespace HDF5

} // namespace dumper
} // namespace akantu
