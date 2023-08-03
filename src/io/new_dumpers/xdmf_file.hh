/**
 * @file   xdmf_file.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Thu Oct 12 2017
 *
 * @brief A Documented file.
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
#include "dumper_file_base.hh"
#include "xml_helper.hh"
/* -------------------------------------------------------------------------- */
#include <filesystem>
#include <sstream>
#include <string_view>
#include <tuple>
#include <typeindex>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

namespace {
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
inline constexpr std::string_view ENDIAN{"Little"};
#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
inline constexpr std::string_view ENDIAN{"Big"};
#else
#error "Does not know from which end to open the eggs"
#endif
} // namespace

#ifndef __AKANTU_XDMF_FILE_HH__
#define __AKANTU_XDMF_FILE_HH__

namespace akantu {
namespace dumper {

  namespace {
    // std::string_view xdmf_type(const ElementType & type) {
    //   static const std::map<ElementType, std::string_view> element_map{
    //       {_point_1, "Polyvertex"},      {_segment_2, "Polyline"},
    //       {_segment_3, "Polyline"},      {_triangle_3, "Triangle"},
    //       {_triangle_6, "Tri_6"},        {_quadrangle_4, "Quadrilateral"},
    //       {_quadrangle_8, "Quad_8"},     {_tetrahedron_4, "Tetrahedron"},
    //       {_tetrahedron_10, "Tet_10"},   {_hexahedron_20, "Hex_20"},
    //       {_hexahedron_8, "Hexahedron"}, {_pentahedron_6, "Wedge"},
    //       {_pentahedron_15, "Wedge_15"}};

    //   return element_map.at(type);
    // }

#define TIDX(x) std::type_index(typeid(x))
    auto xdmf_datatype(const std::type_index & type_idx) {
      static const std::unordered_map<std::type_index,
                                      std::tuple<std::string_view, std::size_t>>
          type_ids{
              {TIDX(int), {"Int", 4}},       {TIDX(unsigned int), {"UInt", 4}},
              {TIDX(int32_t), {"Int", 4}},   {TIDX(int64_t), {"Int", 8}},
              {TIDX(uint32_t), {"UInt", 4}}, {TIDX(uint64_t), {"Int", 8}},
              {TIDX(bool), {"UChar", 1}},    {TIDX(float), {"Float", 4}},
              {TIDX(double), {"Float", 8}},
          };
      return type_ids.at(type_idx);
    }

    std::string xdmf_geometry(char dim) {
      switch (dim) {
      case 1:
        //[[gnu::fallthrougth]];
      case 2:
        return "XY";
      case 3:
        return "XYZ";
      default:
        AKANTU_EXCEPTION("Dim " << dim << " is not a recognized geometry");
      }
    }
  } // namespace

  namespace XDMF {
    namespace fs = std::filesystem;
    class File : public XML::File, public dumper::FileBase {
      using Tag = XML::Tag;
      using CloseTag = XML::CloseTag;

    public:
      /* -------------------------------------------------------------------- */
      File(SupportBase & support, const fs::path & path)
          : XML::File(path), FileBase(support), path(path) {
        *this << "<?xml version=\"1.0\" encoding=\"utf-8\" ?>\n";
        *this << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";

        current_xdmf_element.push_back(1);

        *this << Tag("Xdmf")("Version", "3.0")(
            "xmlns:xi", "http://www.w3.org/2001/XInclude");
        *this << Tag("Domain");
        support.addProperty("xdmf_file", path.string());

        if (not support.hasProperty("dump_count")) {
          return;
        }

        auto mesh_name = support.getName();
        *this << Tag("Grid")("Name", mesh_name + ":all")(
            "GridType", "Collection")("CollectionType", "Temporal");
      }

      /* -------------------------------------------------------------------- */
      File(SupportBase & support, const fs::path & path,
           const std::vector<std::string> & current_xdmf_path,
           const std::vector<int> & current_xdmf_element)
          : XML::File(path), FileBase(support), include(true),
            current_xdmf_path(current_xdmf_path),
            current_xdmf_element(current_xdmf_element), path(path) {}

      /* -------------------------------------------------------------------- */
      XML::File & operator<<(const std::string & str) override {
        return XML::File::operator<<(str);
      }

      /* -------------------------------------------------------------------- */
      XML::File & operator<<(const Tag & tag) override {
        XML::File::operator<<(tag);
        if (tag.autoclosing()) {
          if (not current_xdmf_element.empty()) {
            ++(current_xdmf_element.back());
          }

          return *this;
        }

        auto tag_ = tag.tag();
        if (tag.hasProperty("Name")) {
          tag_ += "[@Name=\"" + tag.getProperty("Name") + "\"]";
        }

        current_xdmf_path.push_back(tag_);
        current_xdmf_element.push_back(1);

        //*this << "<!-- element(" + currentXMFElement() + ") -->\n";
        return *this;
      }

      /* -------------------------------------------------------------------- */
      XML::File & operator<<(const CloseTag & tag) override {
        XML::File::operator<<(tag);
        current_xdmf_path.pop_back();
        current_xdmf_element.pop_back();
        if (not current_xdmf_element.empty()) {
          ++(current_xdmf_element.back());
        }
        return *this;
      }

      /* -------------------------------------------------------------------- */
      void dump() override {
        // if (include) {
        auto count = support.getProperty<Int>("dump_count");
        *this << Tag("Grid")("CollectionType", "Spatial")(
            "GridType", "Collection")("Name", "model_" + std::to_string(count));

        if (support.hasProperty("time")) {
          *this << Tag("Information", true)("Name", "Time")(
              "Value", support.getProperty<Real>("time"));
        }

        FileBase::dump();

        *this << CloseTag{};
        //   return;
        // }

        // auto count = support.getProperty<Int>("dump_count");
        // auto include_path = getIncludeFilename(count);

        // File xdmf_include(support, include_path, current_xdmf_path);
        // xdmf_include.dump();

        // auto sub_path =
        //     fs::relative(include_path, fs::path(path).remove_filename());
        // *this << Tag("xi:include", true)("href", sub_path.string());
      }

    private:
      fs::path getIncludeFilename(std::size_t num) const {
        std::stringstream sstr;
        sstr << std::setw(5) << std::setfill('0') << num;
        auto include_path = path;
        include_path.remove_filename();
        include_path /= "xmfs";
        include_path /=
            path.stem().string() + '_' + sstr.str() + path.extension().string();

        return include_path;
      }

      std::string currentXMFPath() const {
        auto tmp = fs::path("/");
        for (auto && tag : current_xdmf_path) {
          tmp /= tag;
        }
        return tmp.generic_string();
      }

      std::string currentXMFElement() const {
        auto tmp = fs::path("/");
        auto copy = current_xdmf_element;
        copy.pop_back();
        for (auto && tag : copy) {
          tmp /= std::to_string(tag);
        }

        return tmp.generic_string();
      }

      std::string getName(const ID & name, ElementType type,
                          GhostType ghost_type = _not_ghost) const {
        auto full_name = name + ":" + std::to_string(type);
        if (ghost_type != _not_ghost) {
          full_name += ":" + std::to_string(ghost_type);
        }
        return full_name;
      }

    protected:
      template <class T> void xi_include(const T & t, bool parent = false) {
        auto element = t.template getProperty<std::string>("xdmf_element");

        if (parent) {
          element = fs::path(element).parent_path().string();
        }

        *this << Tag("xi:include", true)("xpointer",
                                         "element(" + element + ")");
      }

      /* -------------------------------------------------------------------- */

      /* -------------------------------------------------------------------- */
      void dump(FieldBase & field) override {
        auto && support = field.getSupport();
        auto && tag = Tag("DataItem");

        if (field.hasProperty("attribute_type")) {
          auto && attribute_tag = Tag("Attribute");
          attribute_tag("Name", field.getName());
          attribute_tag("Attribute_Type",
                        field.getProperty<std::string>("attribute_type"));
          attribute_tag("Center",
                        field.getProperty<std::string>("data_location"));

          if (field.hasProperty("finite_element_function")) {
          }

          *this << attribute_tag;
        }

        bool unchanged = true;
        if (field.hasProperty("hdf5_path") and
            field.hasProperty("xdmf_previous_hdf5_path")) {
          unchanged =
              (field.getProperty<std::string>("hdf5_path") ==
               field.getProperty<std::string>("xdmf_previous_hdf5_path"));
        }

        if (field.hasProperty("xdmf_element") and unchanged) {
          xi_include(field);
        } else if (field.hasProperty("hdf5_path")) {
          auto && dims =
              std::to_string(field.getProperty<int>("size")) + " " +
              std::to_string(field.getProperty<int>("nb_components"));
          auto && hdf5_path =
              fs::relative(
                  fs::path(support.getProperty<std::string>("hdf5_file")),
                  fs::path(support.getProperty<std::string>("xdmf_file"))
                      .remove_filename())
                  .string() +
              ":" + field.getProperty<std::string>("hdf5_path");
          auto && [number_type, precision] = xdmf_datatype(field.type());
          /* clang-format off */
          *this << tag("Name", field.getName())
              ("Dimensions", dims)
              ("Format", "HDF")
              ("NumberType", number_type)
              ("Precision", precision)
              ("Endian", ENDIAN);
          /* clang-format on */
          *this << hdf5_path << "\n";

          field.addProperty("xdmf_path", currentXMFPath());
          field.addProperty("xdmf_element", currentXMFElement());
          field.addProperty("xdmf_previous_hdf5_path",
                            field.getProperty<std::string>("hdf5_path"));
          *this << CloseTag{};
        } else {
          AKANTU_EXCEPTION("The XML and Binary format are not implemented yet");
        }

        if (field.hasProperty("attribute_type")) {
          *this << CloseTag{};
        }
      }

      void dump(FieldNodeArrayBase & /*field*/) override {}
      void dump(FieldElementMapArrayBase & /*field*/) override {}

      /* -------------------------------------------------------------------- */
      void dump(Support<Mesh> & support) override {
        auto && mesh_name = support.getName();
        auto && nodes = support.getNodes();

        if (support.hasProperty("xdmf_element") and
            support.getFields().empty()) {
          xi_include(support);
        } else {

          /* clang-format off */
          *this << Tag("Grid")
              ("CollectionType", "Spatial")
              ("GridType", "Collection")
              ("Name", mesh_name);
          /* clang-format on */

          support.addProperty("xdmf_path", currentXMFPath());
          support.addProperty("xdmf_element", currentXMFElement());

          /* clang-format off */
          *this << Tag("Grid")
              ("Name", mesh_name)
              ("GridType", "Uniform");
          /* clang-format on */

          // Topology
          auto & connectivities = support.getConnectivities();
          if (connectivities.hasProperty("xdmf_element")) {
            xi_include(connectivities, true);
          } else {
            UInt nb_elements{0};
            if (connectivities.hasProperty("nb_elements")) {
              nb_elements = connectivities.getProperty<int>("nb_total_element");
            } else {
              nb_elements = connectivities.size().first;
            }

            /* clang-format off */
            *this << Tag("Topology")
                ("TopologyType", "Mixed")
                ("NumberOfElements", nb_elements);
            /* clang-format on */

            dump(aka::as_type<FieldBase>(connectivities));
            *this << CloseTag{};
          }

          // Geometry
          if (nodes.hasProperty("xdmf_element")) {
            xi_include(nodes, true);
          } else {
            *this << Tag("Geometry")(
                "GeometryType", xdmf_geometry(char(nodes.getNbComponent())));
            dump(aka::as_type<FieldBase>(nodes));
            *this << CloseTag{};
          }

          for (auto && data : support.getFields()) {
            auto && field = *data.second;
            switch (field.getFieldType()) {
            case FieldType::_node_array: {
              field.addProperty("attribute_type", "Vector");
              field.addProperty("data_location", "Node");
              this->dump(aka::as_type<FieldBase>(field));
              break;
            }
            case FieldType::_element_map_array: {
              field.addProperty("attribute_type", "Tensor");
              field.addProperty("data_location", "Cell");
              field.addProperty("name", field.getName());
              this->dump(aka::as_type<FieldBase>(field));
              break;
            }
            default:
              AKANTU_EXCEPTION("Unknown field type");
            }
          }

          *this << CloseTag{};
          *this << CloseTag{};
        }

        for (auto && support : support.getSubSupports()) {
          this->dump(*support.second);
        }
      }

      /* -------------------------------------------------------------------- */
      void dump(Support<ElementGroup> & support) override {
        const auto & parent_support = support.getParentSupport();
        fs::path xdmf_parent =
            parent_support.getProperty<std::string>("xdmf_path");
        auto parent_name = parent_support.getName();
        auto name = parent_name + ":" + support.getName();

        if (support.hasProperty("xdmf_element") and
            support.getFields().empty()) {
          xi_include(support);
          return;
        }

        *this << Tag("Grid")("Name", name)("GridType", "Collection")(
            "CollectionType", "Spatial");
        support.addProperty("xdmf_path", currentXMFPath());
        support.addProperty("xdmf_element", currentXMFElement());
        support.addProperty(
            "xdmf_file", parent_support.getProperty<std::string>("xdmf_file"));

        auto & elements = support.getElements();
        for (auto && type : support.elementTypes()) {
          auto xdmf_path = xdmf_parent;
          auto && name = getName(support.getName(), type);

          *this << Tag("Grid")("Name", name)("GridType", "Subset")("Section",
                                                                   "DataItem");
          dump(elements.array(type));
          *this << Tag("Grid")("Name", "Target")("Reference", "XML");
          xdmf_path /= "Grid[@Name=\"" + getName(parent_name, type) + "\"]";
          *this << xdmf_path.generic_string() << "\n";

          *this << CloseTag{};

          for (auto && data : support.getFields()) {
            auto && field = *data.second;
            switch (field.getFieldType()) {
            case FieldType::_node_array: {
              field.addProperty("attribute_type", "Vector");
              field.addProperty("data_location", "Node");
              FileBase::dump(field);
              break;
            }
            case FieldType::_element_map_array: {
              auto && array =
                  aka::as_type<FieldElementMapArrayBase>(field).array(type);
              array.addProperty("attribute_type", "Tensor");
              array.addProperty("data_location", "Cell");
              array.addProperty("name", field.getName());
              dump(array);
              break;
            }
            default:
              AKANTU_EXCEPTION("Unknown field type");
            }
          }

          *this << CloseTag{};
        }

        *this << CloseTag{};
      }

      void dump(SupportBase & support) override { FileBase::dump(support); }

    public:
      bool include{false};
      std::vector<std::string> current_xdmf_path;
      std::vector<int> current_xdmf_element;
      fs::path path;
    };
  } // namespace XDMF
} // namespace dumper
} // namespace akantu

#endif /* __AKANTU_XDMF_FILE_HH__ */
