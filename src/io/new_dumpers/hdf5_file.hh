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
#include "dumper_field.hh"
#include "dumper_file_base.hh"
#include "hdf5_entities.hh"
#include "support.hh"

/* -------------------------------------------------------------------------- */
#include <filesystem>
/* -------------------------------------------------------------------------- */

namespace akantu {
extern std::vector<std::string_view> akantu_dirty_patch;

namespace dumper {

  namespace fs = std::filesystem;
  /* ------------------------------------------------------------------------ */
  class H5File : public FileBase {
  public:
    H5File(SupportBase & support, const fs::path & path);

    ~H5File() override;

    auto filepath() const { return filepath_; }

    bool isAkantuFile();

    void create();
    void open();
    void close();
    /* ---------------------------------------------------------------------- */
    void dump() override;
    void read() override;

  protected:
    template <class T> void writeAttribute(std::string_view name, const T & t) {
      auto & group = *entities.back();
      HDF5::Attribute attr(name, group);
      attr.write(t);
    }

    HDF5::Group & openGroup(const std::string & path,
                            bool create_missing = true);
    HDF5::Symlink & createSymlink(const std::string & path,
                                  const std::string & link);
    HDF5::Dataset & createDataSet(FieldBase & field);
    HDF5::Dataset & openDataSet(FieldBase & field);

  private:
    template <class T> bool unchanged(const T & t) {
      return (t.hasProperty("hdf5_release") and
              t.getRelease() ==
                  t.template getProperty<Release>("hdf5_release"));
    }

    /* ------------------------------------------------------------------ */
    PropertiesManager::slabs_type getSlabs(Support<Mesh> & support);
    PropertiesManager::slabs_type getSlabs(Support<ElementGroup> & support);

    /* -------------------------------------------------------------------- */
    void createDataDescriptions(SupportElements & support);
    /* -------------------------------------------------------------------- */
    void createElementalDataspace(FieldBase & field, HDF5::Dataset & dataset);
    void createNodalDataspace(FieldBase & field, HDF5::Dataset & dataset);
    void createDataspace(FieldBase & field, HDF5::Dataset & dataset);

    /* ---------------------------------------------------------------------- */
    void openElementalDataspace(FieldBase & field, HDF5::Dataset & dataset);
    void openNodalDataspace(FieldBase & field, HDF5::Dataset & dataset);
    void openDataspace(FieldBase & field, HDF5::Dataset & dataset);

    /* ---------------------------------------------------------------------- */
    std::unique_ptr<HDF5::PropertyList> getDataTransferPropertyList();

  protected:
    void dump(FieldNodalArrayBase & field) override;
    void dump(FieldArrayBase & field) override;
    void dump(FieldElementMapArrayBase & field) override;
    void dump(Support<Mesh> & support) override;
    void dump(Support<ElementGroup> & support) override;
    /* ---------------------------------------------------------------------- */
    void read(FieldNodalArrayBase & field) override;
    void read(FieldArrayBase & field) override;
    void read(FieldElementMapArrayBase & field) override;
    void read(Support<Mesh> & support) override;
    void read(Support<ElementGroup> & support) override;

  private:
    std::vector<std::unique_ptr<HDF5::EntityBase>> entities;
    fs::path filepath_;
  };

} // namespace dumper
} // namespace akantu
