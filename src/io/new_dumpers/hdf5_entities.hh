#include "aka_enum_macros.hh"
#if defined(AKANTU_USE_MPI)
#include "mpi_communicator_data.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <filesystem>
#include <hdf5.h>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_HDF5_ENTITIES_HH
#define AKANTU_HDF5_ENTITIES_HH

// #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define ENDIANNESS(x) (x##LE)
// #elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
// #define ENDIANNESS(x) (x##BE)
// #else
// #error "Does not know from which end to open the eggs"
// #endif

namespace akantu {
namespace dumper {
  namespace HDF5 {
    namespace {
      template <class T> std::type_index TIDX(T x) {
        return std::type_index(typeid(x));
      }

      hid_t datatype_id_file(const std::type_index & type_idx) {
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
            {TIDX(double{}), ENDIANNESS(H5T_IEEE_F64)},
            {TIDX(static_cast<const char *>(nullptr)), H5T_C_S1},
            {TIDX(static_cast<char *>(nullptr)), H5T_C_S1},
        };

        if (TIDX(static_cast<const char *>(nullptr)) == type_idx or
            TIDX(static_cast<char *>(nullptr)) == type_idx) {
          auto type_id = H5Tcopy(H5T_C_S1);
          H5Tset_size(type_id, H5T_VARIABLE);
          H5Tset_cset(type_id, H5T_CSET_ASCII);
          return type_id;
        }

        AKANTU_DEBUG_ASSERT(type_ids.find(type_idx) != type_ids.end(),
                            "No HDF5 type known for C++ type "
                                << debug::demangle(type_idx.name()));

        return type_ids.at(type_idx);
      }

      hid_t datatype_id_mem(const std::type_index & type_idx) {
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
            {TIDX(static_cast<const char *>(nullptr)), H5T_C_S1},
            {TIDX(static_cast<char *>(nullptr)), H5T_C_S1},
        };

        if (TIDX(static_cast<const char *>(nullptr)) == type_idx or
            TIDX(static_cast<char *>(nullptr)) == type_idx) {
          auto type_id = H5Tcopy(H5T_C_S1);
          H5Tset_size(type_id, H5T_VARIABLE);
          H5Tset_cset(type_id, H5T_CSET_ASCII);
          return type_id;
        }

        AKANTU_DEBUG_ASSERT(type_ids.find(type_idx) != type_ids.end(),
                            "No HDF5 type known for C++ type "
                                << debug::demangle(type_idx.name()));

        return type_ids.at(type_idx);
      }
    } // namespace

    namespace fs = std::filesystem;
#define AKANTU_ENTITY_TYPES                                                    \
  (group)(dataset)(file)(dataspace)(link)(property)(attribute)

    AKANTU_CLASS_ENUM_DECLARE(EntityType, AKANTU_ENTITY_TYPES)
  } // namespace HDF5
} // namespace dumper

AKANTU_CLASS_ENUM_OUTPUT_STREAM(dumper::HDF5::EntityType, AKANTU_ENTITY_TYPES)
AKANTU_CLASS_ENUM_INPUT_STREAM(dumper::HDF5::EntityType, AKANTU_ENTITY_TYPES)

namespace dumper {
  namespace HDF5 {

    struct EntityBase {
      hid_t id{H5I_INVALID_HID};
      fs::path path;

      EntityBase() = default;
      EntityBase(const fs::path & path) : path(path) {}

      EntityBase(const EntityBase & other) = delete;
      EntityBase(EntityBase && other)
          : id(std::exchange(other.id, H5I_INVALID_HID)),
            path(std::exchange(other.path, {})) {}
      EntityBase & operator=(const EntityBase & other) = delete;
      EntityBase & operator=(EntityBase && other) {
        if (&other != this) {
          id = std::exchange(other.id, H5I_INVALID_HID);
          path = std::exchange(other.path, {});
        }
        return *this;
      }

      virtual ~EntityBase() = default;

      operator hid_t() { return id; }

      virtual void close() {}
    };

    template <EntityType type> struct Entity : public EntityBase {
      using EntityBase::EntityBase;
      using EntityBase::operator=;

      void close() override {
        if (id == H5I_INVALID_HID) {
          return;
        }
        static const std::unordered_map<EntityType, std::function<void(hid_t)>>
            close_func{{EntityType::_group, H5Gclose},
                       {EntityType::_file, H5Fclose},
                       {EntityType::_dataset, H5Dclose},
                       {EntityType::_dataspace, H5Sclose},
                       {EntityType::_attribute, H5Aclose},
                       {EntityType::_property, H5Pclose}};
        if (auto it = close_func.find(type); it != close_func.end()) {
          it->second(id);
        }

        id = H5I_INVALID_HID;
      }

      ~Entity() override { close(); }
    };

    using Symlink = Entity<HDF5::EntityType::_link>;

    /* ---------------------------------------------------------------------- */
    struct PropertyList : public Entity<EntityType::_property> {
      PropertyList(hid_t cls_id) { create(cls_id); }
      PropertyList() { this->id = H5P_DEFAULT; }

      void create(hid_t cls_id) {
        if (id != H5I_INVALID_HID) {
          close();
        }
        id = H5Pcreate(cls_id);
      }

#if defined(AKANTU_USE_MPI)
      void setDxplMPIIO(H5FD_mpio_xfer_t xfer_mode) {
        H5Pset_dxpl_mpio(this->id, xfer_mode);
      }
#endif
    };

    /* ---------------------------------------------------------------------- */
    struct FileAccessPropertyList : public PropertyList {
#if defined(AKANTU_USE_MPI)
      MPI_Info file_info_template{MPI_INFO_NULL};
#endif
      FileAccessPropertyList() : PropertyList(H5P_FILE_ACCESS) {}

#if defined(AKANTU_USE_MPI)
      ~FileAccessPropertyList() override {
        if (file_info_template != MPI_INFO_NULL) {
          MPI_Info_free(&file_info_template);
        }
      }

      void setFaplMPIIO(const Communicator & communicator,
                        const std::map<ID, ID> & infos) {
        if (not infos.empty()) {
          MPI_Info_create(&file_info_template);
          for (auto && [key, value] : infos) {
            MPI_Info_set(file_info_template, key.c_str(), value.c_str());
          }
        }

        auto && mpi_comm = aka::as_type<MPICommunicatorData>(
                               communicator.getCommunicatorData())
                               .getMPICommunicator();

        H5Pset_fapl_mpio(this->id, mpi_comm, file_info_template);
        H5Pset_all_coll_metadata_ops(this->id, true);
        H5Pset_coll_metadata_write(this->id, true);
      }
#endif

      void setLibverBounds(H5F_libver_t low = H5F_LIBVER_LATEST,
                           H5F_libver_t high = H5F_LIBVER_LATEST) {
        H5Pset_libver_bounds(this->id, low, high);
      }
    };

    /* ---------------------------------------------------------------------- */
    struct File : public Entity<EntityType::_file> {
      fs::path filename;
      FileAccessPropertyList access_pl;

      File(const fs::path & filename)
          : Entity<EntityType::_file>("/"), filename(filename) {}

      void open(unsigned mode = H5F_ACC_RDWR) {
        AKANTU_DEBUG_ASSERT(id == H5I_INVALID_HID, "HDF5: File already opened");
        id = H5Fopen(filename.c_str(), mode, access_pl.id);
      }

      void create(unsigned mode = H5F_ACC_RDWR, hid_t fcpl_id = H5P_DEFAULT) {
        AKANTU_DEBUG_ASSERT(id == H5I_INVALID_HID, "HDF5: File already opened");
        id = H5Fcreate(filename.c_str(), mode, fcpl_id, access_pl.id);
      }

      void close() override {
#if !defined(AKANTU_NDEBUG)
        if (id == H5I_INVALID_HID) {
          return;
        }

        auto nids = H5Fget_obj_count(this->id, H5F_OBJ_ALL);

        std::vector<hid_t> open_ids(nids);
        H5Fget_obj_ids(this->id, H5F_OBJ_ALL, nids, open_ids.data());

        const std::unordered_map<H5I_type_t, ID> names{
            {H5I_DATATYPE, "datatype"}, {H5I_DATASPACE, "dataspace"},
            {H5I_DATASET, "dataset"},   {H5I_GROUP, "group"},
            {H5I_FILE, "file"},         {H5I_ATTR, "attribute"},
            {H5I_BADID, "bad_id"}};

        for (auto && oid : open_ids) {
          auto type = H5Iget_type(oid);
          if (type != H5I_FILE) {
            AKANTU_DEBUG_WARNING("HDF5 forgot to close id: "
                                 << oid << " of type: " << names.at(type));
          }
        }
#endif
        Entity<EntityType::_file>::close();
      }

      auto & getAccessPropertyList() { return access_pl; }
    };

    /* ---------------------------------------------------------------------- */
    struct Group;

    template <EntityType type> struct GenericGroup : public Entity<type> {
      const EntityBase & parent;
      static_assert(type == EntityType::_group or type == EntityType::_dataset);

      GenericGroup(const std::string & path, const EntityBase & parent);
      void open(hid_t apl_id = H5P_DEFAULT);
      bool exists(hid_t lapl_id = H5P_DEFAULT);

      virtual void create(hid_t lcpl_id = H5P_DEFAULT,
                          hid_t gcpl_id = H5P_DEFAULT,
                          hid_t gapl_id = H5P_DEFAULT) = 0;

      void createOrOpen() {
        if (not this->exists()) {
          AKANTU_DEBUG_INFO("HDF5: Creating " << type << " "
                                              << this->path.generic_string());
          this->create();
        } else {
          AKANTU_DEBUG_INFO("HDF5: Opening existing "
                            << type << " " << this->path.generic_string());
          this->open();
        }
      }
    };

    /* ---------------------------------------------------------------------- */
    struct Group : public GenericGroup<EntityType::_group> {
      using GenericGroup<EntityType::_group>::GenericGroup;

      void create(hid_t lcpl_id = H5P_DEFAULT, hid_t gcpl_id = H5P_DEFAULT,
                  hid_t gapl_id = H5P_DEFAULT) override {
        AKANTU_DEBUG_ASSERT(this->id == H5I_INVALID_HID,
                            "HDF5: Group already opened");
        this->id =
            H5Gcreate(parent.id, this->path.c_str(), lcpl_id, gcpl_id, gapl_id);
      }
    };

    /* ---------------------------------------------------------------------- */
    template <EntityType type>
    GenericGroup<type>::GenericGroup(const std::string & path,
                                     const EntityBase & parent)
        : Entity<type>(parent.path / path), parent(parent) {}

    template <EntityType type> void GenericGroup<type>::open(hid_t apl_id) {
      AKANTU_DEBUG_ASSERT(this->id == H5I_INVALID_HID,
                          "HDF5: GenericGroup(" << type << ") already opened");

      if (not exists()) {
        AKANTU_EXCEPTION("HDF5: " << type << ": " << this->path
                                  << " does not exists");
      }

      if constexpr (type == EntityType::_group) {
        this->id = H5Gopen(this->parent.id, this->path.c_str(), apl_id);
      } else if constexpr (type == EntityType::_dataset) {
        this->id = H5Dopen(this->parent.id, this->path.c_str(), apl_id);
      }
    }

    template <EntityType type> bool GenericGroup<type>::exists(hid_t lapl_id) {
      auto status = H5Lexists(this->parent.id, this->path.c_str(), lapl_id);
      return (status > 0);
    }

    /* ---------------------------------------------------------------------- */
    struct Dataspace : public Entity<EntityType::_dataspace> {
      Dataspace(const std::array<hsize_t, 2> & dims) {
        id = H5Screate_simple(dims.size(), dims.data(), nullptr);
      }

      Dataspace(const Dataspace & other) { id = H5Scopy(other.id); }
    };

    /* ---------------------------------------------------------------------- */
    struct Dataset : public GenericGroup<EntityType::_dataset> {
      std::unique_ptr<Dataspace> dataspace;
      std::unique_ptr<Dataspace> memoryspace;
      std::unique_ptr<Dataspace> filespace;
      std::type_index datatype;

      Dataset(const std::string & path, const Group & parent,
              const std::type_index & datatype)
          : GenericGroup<EntityType::_dataset>(path, parent),
            datatype(datatype) {}

      void create(hid_t lcpl_id = H5P_DEFAULT, hid_t dcpl_id = H5P_DEFAULT,
                  hid_t dapl_id = H5P_DEFAULT) override {
        AKANTU_DEBUG_ASSERT(this->id == H5I_INVALID_HID,
                            "HDF5: Dataset already opened");
        this->id = H5Dcreate(this->parent.id, this->path.c_str(),
                             datatype_id_mem(datatype), this->dataspace->id,
                             lcpl_id, dcpl_id, dapl_id);
      }

      void write(const void * data, hid_t dxpl_id = H5P_DEFAULT) {
        AKANTU_DEBUG_INFO("HDF5: Writing dataset " << this->path);
        H5Dwrite(this->id, datatype_id_file(datatype), memoryspace->id,
                 filespace->id, dxpl_id, data);
      }

      void read(void * data, hid_t dxpl_id = H5P_DEFAULT) {
        AKANTU_DEBUG_INFO("HDF5: Writing dataset " << this->path);
        H5Dread(this->id, datatype_id_mem(datatype), memoryspace->id,
                filespace->id, dxpl_id, data);
      }
    };

    /* ---------------------------------------------------------------------- */
    struct Attribute : public Entity<EntityType::_attribute> {
      const EntityBase & location;
      ID name;
      Attribute(std::string_view name, const EntityBase & location)
          : Entity<EntityType::_attribute>(location.path), location(location),
            name(name) {}

      [[nodiscard]] bool exists() const {
        return (H5Aexists(location.id, name.data()) > 0);
      }

    private:
      template <class T>
      void read(T * data, hid_t type_id_mem = H5I_INVALID_HID) {
        if (this->id == H5I_INVALID_HID) {
          open();
        }

        AKANTU_DEBUG_INFO("HDF5: Reading attribute " << (this->path / name));

        auto type_id = type_id_mem;
        if (type_id == H5I_INVALID_HID) {
          type_id = datatype_id_mem(TIDX(T{}));
        }

        H5Aread(this->id, type_id, data);

        // if (type_id_mem != type_id) {
        //   H5Tclose(type_id);
        // }
      }

      template <typename T, std::enable_if_t<std::is_scalar_v<T>> * = nullptr>
      void write(T * data, hsize_t dim, hid_t acpl_id = H5P_DEFAULT,
                 hid_t aapl_id = H5P_DEFAULT) {
        auto type_id = datatype_id_file(TIDX(T{}));
        auto space_id = H5Screate_simple(1, &dim, &dim);

        AKANTU_DEBUG_INFO("HDF5: Writing attribute " << (this->path / name));

        create(type_id, space_id, acpl_id, aapl_id);

        H5Awrite(this->id, type_id, data);
      }

    public:
      template <typename T, std::enable_if_t<std::is_scalar_v<T>> * = nullptr>
      void write(const T & t, hid_t acpl_id = H5P_DEFAULT,
                 hid_t aapl_id = H5P_DEFAULT) {
        write(&t, 1, acpl_id, aapl_id);
      }

      void write(std::string_view value, hid_t acpl_id = H5P_DEFAULT,
                 hid_t aapl_id = H5P_DEFAULT) {
        const char * str = value.data();
        write(&str, 1, acpl_id, aapl_id);
      }

      template <typename T, std::enable_if_t<std::is_scalar_v<T>> * = nullptr>
      void write(const std::vector<T> & t, hid_t acpl_id = H5P_DEFAULT,
                 hid_t aapl_id = H5P_DEFAULT) {
        write(t.data(), t.size(), acpl_id, aapl_id);
      }

      void write(const std::vector<std::string_view> & ts,
                 hid_t acpl_id = H5P_DEFAULT, hid_t aapl_id = H5P_DEFAULT) {
        std::vector<const char *> strs(ts.size());
        for (auto && [str, t] : zip(strs, ts)) {
          str = t.data();
        }
        write(strs, acpl_id, aapl_id);
      }

      template <class T, std::enable_if_t<std::is_scalar_v<T>> * = nullptr>
      T read() {
        T t;
        read(&t);
        return t;
      }

      template <class T> void read(std::vector<T> & ts) {
        open();
        auto space_id = H5Aget_space(this->id);
        auto ndims = H5Sget_simple_extent_ndims(space_id);
        if (ndims != 1) {
          AKANTU_EXCEPTION("Cannot read the attribute "
                           << location.path << "/" << name
                           << " the dimensions do not match");
        }

        hsize_t dim, maxdim;
        H5Sget_simple_extent_dims(space_id, &dim, &maxdim);
        ts.resize(dim);

        if constexpr (std::is_same_v<T, std::string>) {
          std::vector<char *> strs(ts.size());
          auto type_id_mem = datatype_id_mem(TIDX(strs[0]));

          read(strs.data(), type_id_mem);

          for (auto && [str, t] : zip(strs, ts)) {
            t = std::string(str);
          }
          H5Dvlen_reclaim(type_id_mem, space_id, H5P_DEFAULT, strs.data());
          H5Tclose(type_id_mem);
        } else {
          read(ts.data());
        }
      }

      template <class T,
                std::enable_if_t<std::is_same_v<T, std::string>> * = nullptr>
      T read() {
        std::vector<std::string> ts;
        read(ts);
        return ts[0];
      }

    private:
      auto open(hid_t aapl_id = H5P_DEFAULT) {
        if (not exists()) {
          AKANTU_EXCEPTION("HDF5: Attribute: "
                           << name << " does not exists at location "
                           << location.path.generic_string());
        }

        this->id = H5Aopen(location.id, name.data(), aapl_id);
      }

      void create(hid_t type_id, hid_t space_id, hid_t acpl_id = H5P_DEFAULT,
                  hid_t aapl_id = H5P_DEFAULT) {
        if (H5Aexists(location.id, name.c_str()) > 0) {
          H5Adelete(location.id, name.c_str());
        }

        this->id = H5Acreate(location.id, name.c_str(), type_id, space_id,
                             acpl_id, aapl_id);
      }
    };

    extern "C" {
    inline herr_t hdf5_error_walk(unsigned int n, const H5E_error2_t * err_desc,
                                  void * client_data) {
      const int MSG_SIZE = 64;
      auto & messages =
          *reinterpret_cast<std::vector<std::string> *>(client_data);
      char maj[MSG_SIZE];
      char min[MSG_SIZE];
      char cls[MSG_SIZE];
      char file[MSG_SIZE];
      char func[MSG_SIZE];
      char desc[MSG_SIZE];

      std::strncpy(file, err_desc->file_name, MSG_SIZE);
      std::strncpy(func, err_desc->func_name, MSG_SIZE);
      std::strncpy(desc, err_desc->desc, MSG_SIZE);

      if (H5Eget_class_name(err_desc->cls_id, cls, MSG_SIZE) < 0) {
        return -1;
      }

      if (H5Eget_msg(err_desc->maj_num, NULL, maj, MSG_SIZE) < 0) {
        return -1;
      }

      if (H5Eget_msg(err_desc->min_num, NULL, min, MSG_SIZE) < 0) {
        return -1;
      }

      char stream[7 * MSG_SIZE];
      std::snprintf(stream, 7 * MSG_SIZE,
                    "%s [error: #%03d] - %s (%s: %s) in %s(): [%s:%u]", cls, n,
                    maj, min, desc, func, file, err_desc->line);
      messages.emplace_back(stream);
      return 0;
    }

    inline herr_t hdf5_error_handler(hid_t estack_id, void * client_data) {
      auto & file = *reinterpret_cast<fs::path *>(client_data);

      std::vector<std::string> messages;
      H5Ewalk(estack_id, H5E_WALK_DOWNWARD, hdf5_error_walk, &messages);

      std::stringstream sstr;
      sstr << "HDF5 errors while accessing file \"" << file.generic_string()
           << "\":";
      for (auto && msg : messages) {
        sstr << "\n" << std::string(4, AKANTU_INDENT) << msg;
      }

      AKANTU_EXCEPTION(sstr.str());

      return 0;
    }
    }
  } // namespace HDF5
} // namespace dumper
} // namespace akantu

#endif // AKANTU_HDF5_ENTITIES_HH
