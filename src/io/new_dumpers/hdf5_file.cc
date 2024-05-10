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
#include "dumper_types.hh"
/* -------------------------------------------------------------------------- */
#include "aka_iterators.hh"
#include "dumper_field.hh"
#include "dumper_file_base.hh"
#include "hdf5_entities.hh"
#include "hdf5_file.hh"
#include "support.hh"
/* -------------------------------------------------------------------------- */
#include <array>
#include <filesystem>
#include <sstream>
#include <typeindex>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

namespace akantu {
extern std::vector<std::string_view> akantu_dirty_patch;

namespace dumper {
  namespace fs = std::filesystem;

  namespace HDF5 {
    /* ---------------------------------------------------------------------- */
    template <class T>
    class HDF5NodeFieldTemplateIterator
        : public FieldComputeArray<T, T, FieldNodalArrayBase> {
      using parent = FieldComputeArray<T, T, FieldNodalArrayBase>;
      using in_field = FieldArrayTemplateBase<T, FieldNodalArrayBase>;

    public:
      using parent::parent;

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
            make_view(const_cast<const in_field &>(*this->array_in).getArray(),
                      nb_component)
                .begin();
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

      void back_compute() override {
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
            *in_it = *out_it;
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
  } // namespace HDF5

  /* ------------------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */
  H5File::H5File(SupportBase & support, const fs::path & path)
      : FileBase(support), filepath_(path) {
    auto path_wof = path;
    path_wof.remove_filename();
    if (not fs::exists(path_wof)) {
      fs::create_directories(path_wof);
    }

    H5Eset_auto(H5E_DEFAULT, HDF5::hdf5_error_handler, &filepath_);

    auto file = std::make_unique<HDF5::File>(path);

    auto & fapl = file->getAccessPropertyList();
#if defined(AKANTU_USE_MPI)
    fapl.setFaplMPIIO(support.getCommunicator(),
                      {{"access_style", "write_once"},
                       {"collective_buffering", "true"},
                       {"cb_block_size", "1048576"},
                       {"cb_buffer_size", "4194304"}});
#endif

    fapl.setLibverBounds();

    entities.push_back(std::move(file));
  }

  /* ------------------------------------------------------------------------ */
  H5File::~H5File() { close(); }

  /* ------------------------------------------------------------------------ */
  void H5File::close() {
    AKANTU_DEBUG_ASSERT(entities.size() == 1,
                        "Not all object where closed, or too many");
    entities[0]->close();
  }

  /* ------------------------------------------------------------------------ */
  void H5File::create() {
    auto & file = aka::as_type<HDF5::File>(*entities[0]);
    file.create(H5F_ACC_TRUNC, H5P_DEFAULT);

    support.addProperty("hdf5_file", filepath().string());

    openGroup("metadata");
    writeAttribute("library", "akantu");
    writeAttribute("version_akantu", getVersion());
    if (not akantu_dirty_patch.empty()) {
      writeAttribute("akantu_patch", akantu_dirty_patch);
    }
    writeAttribute("version_file", 1);
    writeAttribute("seed", debug::_global_seed);
    writeAttribute("nb_proc", support.getCommunicator().getNbProc());
    entities.pop_back();
  }

  /* ------------------------------------------------------------------------ */
  void H5File::open() {
    auto & file = aka::as_type<HDF5::File>(*entities[0]);
    file.open(H5F_ACC_RDWR);
  }

  /* ------------------------------------------------------------------------ */
  bool H5File::isAkantuFile() {
    auto && parent = *entities.back();
    auto && group = std::make_unique<HDF5::Group>("metadata", parent);
    if (not group->exists()) {
      return false;
    }

    group->open();
    auto && attribute = std::make_unique<HDF5::Attribute>("library", *group);
    if (not attribute->exists()) {
      return false;
    }

    auto value = attribute->read<std::string>();
    return value == "akantu";
  }

  /* ------------------------------------------------------------------------ */
  HDF5::Group & H5File::openGroup(const std::string & path,
                                  bool create_missing) {
    auto & parent = *entities.back();

    auto && new_group = std::make_unique<HDF5::Group>(path, parent);
    if (create_missing) {
      new_group->createOrOpen();
    } else {
      new_group->open();
    }

    entities.push_back(std::move(new_group));
    return aka::as_type<HDF5::Group>(*entities.back());
  }

  /* ------------------------------------------------------------------------ */
  HDF5::Symlink & H5File::createSymlink(const std::string & path,
                                        const std::string & link) {
    auto & group = *entities.back();

    auto && new_link = std::make_unique<HDF5::Symlink>(group.path / path);

    auto status [[maybe_unused]] =
        H5Lcreate_soft(link.c_str(), group.id, new_link->path.c_str(),
                       H5P_DEFAULT, H5P_DEFAULT);

    AKANTU_DEBUG_ASSERT(status >= 0, "Could not create HDF5 link "
                                         << new_link->path.c_str() << " -> "
                                         << link.c_str()
                                         << ", status=" << status);

    entities.push_back(std::move(new_link));
    return aka::as_type<HDF5::Symlink>(*entities.back());
  }

  /* ------------------------------------------------------------------------ */
  HDF5::Dataset & H5File::createDataSet(FieldBase & field) {
    auto & group = aka::as_type<HDF5::Group>(*entities.back());
    auto data_set =
        std::make_unique<HDF5::Dataset>(field.getName(), group, field.type());

    createDataspace(field, *data_set);

    data_set->createOrOpen();

    field.addProperty("hdf5_path", data_set->path.generic_string());

    entities.push_back(std::move(data_set));
    return aka::as_type<HDF5::Dataset>(*entities.back());
  }

  /* ------------------------------------------------------------------ */
  PropertiesManager::slabs_type H5File::getSlabs(Support<Mesh> & support) {
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
  PropertiesManager::slabs_type
  H5File::getSlabs(Support<ElementGroup> & support) {
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

  /* ------------------------------------------------------------------------ */
  void H5File::createDataDescriptions(SupportElements & support) {
    auto & support_base = dynamic_cast<SupportBase &>(support);

    // Information needed for elements
    if (not support_base.hasProperty("hdf5_description_release") or
        support_base.getProperty<Release>("hdf5_description_release") !=
            support_base.getRelease()) {
      support.updateOffsets();
    }

    auto slabs = tuple_dispatch<AllSupportTypes>(
        [&](auto type) { return getSlabs(support_cast(support_base, type)); },
        support_base.getType());

    support_base.addProperty("hdf5_nodes_slabs", slabs);
    support_base.addProperty("hdf5_description_release",
                             support_base.getRelease());
  }

  /* ------------------------------------------------------------------------ */
  void H5File::createElementalDataspace(FieldBase & field,
                                        HDF5::Dataset & dataset) {
    auto && field_element = aka::as_type<FieldElementalArrayBase>(field);

    auto && support = dynamic_cast<const SupportElements &>(field.getSupport());
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

  /* ------------------------------------------------------------------------ */
  void H5File::createNodalDataspace(FieldBase & field,
                                    HDF5::Dataset & dataset) {
    auto && support = dynamic_cast<const SupportElements &>(field.getSupport());

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

  /* ---------------------------------------------------------------------- */
  void H5File::createDataspace(FieldBase & field, HDF5::Dataset & dataset) {
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

  /* ------------------------------------------------------------------------ */
  std::unique_ptr<HDF5::PropertyList> H5File::getDataTransferPropertyList() {
    auto dxpl = std::make_unique<HDF5::PropertyList>();
#if defined(AKANTU_USE_MPI)
    dxpl->create(H5P_DATASET_XFER);
    //        H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
    dxpl->setDxplMPIIO(H5FD_MPIO_INDEPENDENT);
#endif
    return dxpl;
  }

  /* ------------------------------------------------------------------------ */
  void H5File::dump(FieldNodalArrayBase & field) {
    auto hdf5_field = HDF5::make_hdf5_field(field);
    dump(dynamic_cast<FieldArrayBase &>(*hdf5_field));
  }

  /* ------------------------------------------------------------------------ */
  void H5File::dump(FieldArrayBase & field) {
    field.addProperty("size", field.size());
    field.addProperty("nb_components", field.getNbComponent());

    auto & dataset = createDataSet(field);
    auto dxpl = getDataTransferPropertyList();

    field.update();

    auto && cfield = const_cast<const FieldArrayBase &>(field);

    dataset.write(cfield.data(), *dxpl);

    field.addProperty("hdf5_release", cfield.getRelease());

    entities.pop_back(); // dataset
  }

  /* ------------------------------------------------------------------------ */
  void H5File::dump(FieldElementMapArrayBase & field) {
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
        array.addProperty("name", field.getName() + ":" + std::to_string(type));
      }
      dump(array);
    }
    field.addProperty(
        "hdf5_release",
        const_cast<const FieldElementMapArrayBase &>(field).getRelease());
    entities.pop_back();
  }

  /* ------------------------------------------------------------------------ */
  void H5File::dump(Support<Mesh> & support) {
    openGroup(support.getName());

    if (unchanged(support)) {
      createSymlink("topology", support.getProperty<std::string>("hdf5_path"));
      entities.pop_back();
    } else {
      createDataDescriptions(support);

      auto & topology = openGroup("topology");
      support.addProperty("hdf5_path", topology.path.generic_string());

      openGroup("nodes");
      auto && nodes = support.getNodes();
      dump(nodes);
      entities.pop_back();

      auto && connectivities =
          make_field(support.getConnectivities().getSharedPointer(), support,
                     toVTKConnectivity(support.getMesh()),
                     toAkantuConnectivity(support.getMesh()));
      dump(*connectivities);
      entities.pop_back();

      support.addProperty("hdf5_release", support.getRelease());
    }

    openGroup("groups");
    for (auto && [_, support] : support.getSubSupports()) {
      createDataDescriptions(aka::as_type<SupportElements>(*support));
      FileBase::dump(*support);
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

  /* ------------------------------------------------------------------------ */
  void H5File::dump(Support<ElementGroup> & support) {
    const auto & parent_support = support.getParentSupport();
    openGroup(support.getName());

    if (unchanged(support)) {
      createSymlink("topology", support.getProperty<std::string>("hdf5_path"));
      entities.pop_back();
    } else {
      support.addProperty("hdf5_file",
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

  /* ------------------------------------------------------------------------ */
  void H5File::dump() {
    openGroup("steps");

    if (support.hasProperty("time")) {
      openGroup(std::to_string(support.getProperty<Int>("dump_count")));
      writeAttribute("time", support.getProperty<double>("time"));
    }

    FileBase::dump(support);

    if (support.hasProperty("time")) {
      entities.pop_back();
    }
    entities.pop_back();

    if (support.hasProperty("time")) {
      openGroup("metadata");
      writeAttribute("last_step", support.getProperty<Int>("dump_count"));
      entities.pop_back();
    }
  }

  /* ------------------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */
  void H5File::openElementalDataspace(FieldBase & field,
                                      HDF5::Dataset & dataset) {
    auto && field_element = aka::as_type<FieldElementalArrayBase>(field);

    auto && support = dynamic_cast<const SupportElements &>(field.getSupport());
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

  /* ------------------------------------------------------------------------ */
  void H5File::openNodalDataspace(FieldBase & field, HDF5::Dataset & dataset) {
    auto && support = dynamic_cast<const SupportElements &>(field.getSupport());

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

  /* ------------------------------------------------------------------------ */
  void H5File::openDataspace(FieldBase & field, HDF5::Dataset & dataset) {
    field.addProperty<hid_t>("hdf5_filespace", H5I_INVALID_HID);

    using dumper::FieldType;
    switch (field.getFieldType()) {
    case FieldType::_node_array:
      openNodalDataspace(field, dataset);
      break;
    case FieldType::_element_array:
      openElementalDataspace(field, dataset);
      break;
    default:
      AKANTU_EXCEPTION("The field type is not properly defined");
      break;
    }
  }

  /* ---------------------------------------------------------------------- */
  HDF5::Dataset & H5File::openDataSet(FieldBase & field) {
    auto & group = aka::as_type<HDF5::Group>(*entities.back());
    auto data_set =
        std::make_unique<HDF5::Dataset>(field.getName(), group, field.type());

    data_set->open();
    openDataspace(field, *data_set);

    entities.push_back(std::move(data_set));
    return aka::as_type<HDF5::Dataset>(*entities.back());
  }

  /* ------------------------------------------------------------------------ */
  void H5File::read() {
    auto & group = openGroup("metadata", false);
    auto && attribute = std::make_unique<HDF5::Attribute>("last_step", group);
    auto last_step_exists = attribute->exists();
    Int last_step{-1};
    if (last_step_exists) {
      last_step = attribute->read<Int>();
    }
    entities.pop_back(); // metadata group

    openGroup("steps", false);

    if (last_step_exists) {
      openGroup(std::to_string(last_step), false);
    }

    FileBase::read(support);

    if (last_step_exists) {
      entities.pop_back(); // last_step group
    }
    entities.pop_back(); // steps group
  }

  /* ------------------------------------------------------------------------ */
  void H5File::read(Support<Mesh> & support) {
    openGroup(support.getName(), false);

    // Do not read the mesh topology, it is supposed to already exists

    // readDataDescriptions(support);

    openGroup("groups", false);
    for (auto && [_, support] : support.getSubSupports()) {
      // readDataDescriptions(aka::as_type<SupportElements>(*support));
      FileBase::read(*support);
    }
    // close groups
    entities.pop_back();

    openGroup("data", false);
    for (auto && [_, field] : support.getFields()) {
      if (not field->is_restart()) {
        FileBase::dump(*field);
      }
    }
    entities.pop_back();

    // close mesh
    entities.pop_back();
  }

  /* ------------------------------------------------------------------------ */
  void H5File::read(Support<ElementGroup> & support) {
    openGroup(support.getName(), false);

    // Do not read the group topology, it is supposed to already exists.

    openGroup("data", false);
    for (auto && [_, field] : support.getFields()) {
      if (not field->is_restart()) {
        FileBase::dump(*field);
      }
    }

    entities.pop_back();
    entities.pop_back();
  }

  /* ------------------------------------------------------------------------ */
  void H5File::read(FieldNodalArrayBase & field) {
    auto hdf5_field = HDF5::make_hdf5_field(field);
    read(dynamic_cast<FieldArrayBase &>(*hdf5_field));
  }

  /* ------------------------------------------------------------------------ */
  void H5File::read(FieldArrayBase & field) {
    field.addProperty("size", field.size());
    field.addProperty("nb_components", field.getNbComponent());

    auto & dataset = openDataSet(field);
    auto dxpl = getDataTransferPropertyList();

    dataset.read(field.data(), *dxpl);
    field.back_propagate();

    entities.pop_back(); // dataset
  }

  /* ------------------------------------------------------------------------ */
  void H5File::read(FieldElementMapArrayBase & field) {
    openGroup(field.getName(), false);
    for (auto && type : field.elementTypes(_not_ghost)) {
      auto & array = field.array(type);
      read(array);
    }
    entities.pop_back();
  }

} // namespace dumper
} // namespace akantu
