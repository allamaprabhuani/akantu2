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
    std::string_view xdmf_type(const ElementType & type) {
      static const std::map<ElementType, std::string_view> element_map{
          {_point_1, "Polyvertex"},      {_segment_2, "Polyline"},
          {_segment_3, "Polyline"},      {_triangle_3, "Triangle"},
          {_triangle_6, "Tri_6"},        {_quadrangle_4, "Quadrilateral"},
          {_quadrangle_8, "Quad_8"},     {_tetrahedron_4, "Tetrahedron"},
          {_tetrahedron_10, "Tet_10"},   {_hexahedron_20, "Hex_20"},
          {_hexahedron_8, "Hexahedron"}, {_pentahedron_6, "Wedge"},
          {_pentahedron_15, "Wedge_15"}};

      return element_map.at(type);
    }

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

    std::tuple<Int, std::string>
    xdmf_finite_element_function(ElementType type) {
      std::unordered_map<ElementType, std::tuple<Int, std::string>> infos{
          {_segment_2, {1, "interval"}},
          {_segment_3, {2, "interval"}},
          {_triangle_3, {1, "triangle"}},
          {_triangle_6, {2, "triangle"}},
          {_quadrangle_4, {1, "quadrilateral"}},
          {_quadrangle_8, {2, "quadrilateral"}},
          {_tetrahedron_4, {1, "tetrahedron"}},
          {_tetrahedron_10, {2, "tetrahedron"}},
          {_hexahedron_20, {1, "hexahedron"}},
          {_hexahedron_8, {2, "hexahedron"}}};
      return infos.at(type);
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
        auto count = support.getProperty<Int>("dump_count");
        *this << Tag("Grid")("CollectionType", "Spatial")(
            "GridType", "Collection")("Name", "model_" + std::to_string(count));

        if (support.hasProperty("time")) {
          *this << Tag("Information", true)("Name", "Time")(
              "Value", support.getProperty<Real>("time"));
        }

        FileBase::dump();

        *this << CloseTag{};
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

      template <class T> bool unchanged(const T & t) {
        if (t.hasProperty("hdf5_path") and
            t.hasProperty("xdmf_previous_hdf5_path")) {
          return (
              t.template getProperty<std::string>("hdf5_path") ==
              t.template getProperty<std::string>("xdmf_previous_hdf5_path"));
        }
        return false;
      }

      /* -------------------------------------------------------------------- */
      [[nodiscard]] std::string
      attributeLocation(const FieldBase & field) const {
        switch (field.getFieldType()) {
        case FieldType::_node_array_function: /* FALLTHRU */
        case FieldType::_node_array: {
          return "Node";
        }
        case FieldType::_element_map_array_function: /* FALLTHRU */
        case FieldType::_element_map_array: {
          return "Cell";
        }
        case FieldType::_internal_field_function: /* FALLTHRU */
        case FieldType::_internal_field: {
          return "Cell";
        }
        default:
          AKANTU_EXCEPTION("Unknown field type");
        }
      }

      /* -------------------------------------------------------------------- */
      void dumpAttribute(FieldBase & field) {
        if (not field.hasProperty("data_location")) {
          field.addProperty("data_location", attributeLocation(field));
        }

        if (field.hasProperty("xdmf_element") and unchanged(field)) {
          xi_include(field, true);
        } else {
          auto && attribute_tag = Tag("Attribute");
          attribute_tag("Name", field.getName());
          if (field.hasProperty("attribute_type")) {
            attribute_tag("Attribute_Type",
                          field.getProperty<std::string>("attribute_type"));
          }
          attribute_tag("Center",
                        field.getProperty<std::string>("data_location"));
          if (field.hasProperty("finite_element_function")) {
            auto [type, ghost_type] =
                field.getProperty<std::tuple<ElementType, GhostType>>(
                    "element_type");
            auto [degree, family] = xdmf_finite_element_function(type);

            // clang-format off
            attribute_tag("ItemType", "FiniteElementFunction")
                ("ElementFamily", "DG")
                ("ElementDegree", degree)
                ("ElementCell", family);
            // clang-format on
          }

          *this << attribute_tag;

          if (field.hasProperty("finite_element_function")) {
            auto [type, ghost_type] =
                field.getProperty<std::tuple<ElementType, GhostType>>(
                    "element_type");
            auto & xdmf_connectivities =
                aka::as_type<FieldElementMapArrayBase>(
                    *field.getSupport().getField("xdmf_connectivities"))
                    .array(type, ghost_type);
            if (xdmf_connectivities.hasProperty("xdmf_element") and
                unchanged(xdmf_connectivities)) {
              xi_include(xdmf_connectivities);
            } else {
              dump(xdmf_connectivities);
            }
          }

          this->dump(field);

          *this << CloseTag{};
        }
      }

      void dumpAttribute(FieldBase & field_el, ElementType type,
                         GhostType ghost_type = _not_ghost) {
        auto & field = aka::as_type<FieldElementMapArrayBase>(field_el).array(
            type, ghost_type);

        if (field.hasProperty("xdmf_element") and unchanged(field)) {
          xi_include(field, true);
        } else {
          auto name = field.getName();

          // if (is_quadrature_points_field(field_el)) {
          //   field.addProperty("finite_element_function", true);
          // }

          field.addProperty("element_type", std::pair{type, ghost_type});
          field.addProperty("name", field_el.getName());
          field.addProperty("data_location", attributeLocation(field_el));
          dumpAttribute(field);
          field.addProperty("name", name);
        }
      }

      /* -------------------------------------------------------------------- */
      void dump(FieldBase & field) override { FileBase::dump(field); }

      /* -------------------------------------------------------------------- */
      void dump(FieldArrayBase & field) override {
        auto && support = field.getSupport();
        auto && tag = Tag("DataItem");

        if (field.hasProperty("xdmf_element") and unchanged(field) and
            not field.hasProperty("xdmf_no_link")) {
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
      }

      void dump(FieldElementMapArrayBase & /*field*/) override {}

      void dumpTopology(FieldArrayBase & connectivities, ElementType type) {
        if (connectivities.hasProperty("xdmf_element")) {
          xi_include(connectivities, true);
        } else {
          auto nb_elements = connectivities.size();
          /* clang-format off */
          *this << Tag("Topology")
              ("TopologyType", xdmf_type(type))
              ("NumberOfElements", nb_elements);
          /* clang-format on */

          dump(connectivities);
          *this << CloseTag{};
        }
      }

      void dumpGeometry(FieldArrayBase & nodes) {
        if (nodes.hasProperty("xdmf_element")) {
          xi_include(nodes, true);
        } else {
          auto geom_type = xdmf_geometry(char(nodes.getNbComponent()));
          *this << Tag("Geometry")("GeometryType", geom_type);
          dump(nodes);
          *this << CloseTag{};
        }
      }

      void dumpSet(Support<ElementGroup> & support, ElementType type,
                   GhostType ghost_type = _not_ghost) {
        const auto & parent_support = support.getParentSupport();
        support.addProperty(
            "xdmf_file", parent_support.getProperty<std::string>("xdmf_file"));

        if (not support.contains(type, ghost_type)) {
          return;
        }

        Int nb_fields{};
        bool fields_unchanged = true;
        for (auto && [_, field] : support.getFields()) {
          if (field->is_elemental_field()) {
            nb_fields++;
            fields_unchanged &= unchanged(*field);
          }
        }
        if (nb_fields == 0) {
          return;
        }

        auto && elements = support.getElements().array(type, ghost_type);

        if (elements.hasProperty("xdmf_element") and unchanged(elements) and
            fields_unchanged) {
          xi_include(elements, true);
          return;
        }

        /* clang-format off */
        *this << Tag("Set")
            ("Name", support.getName() + ":" + std::to_string(type))
            ("SetType", "Cell");
        /* clang-format on */
        elements.addProperty("xdmf_no_link", true);
        dump(elements);

        for (auto && [_, field] : support.getFields()) {
          if (not field->is_elemental_field()) {
            continue;
          }

          if (not field->is_visualization()) {
            continue;
          }
          dumpAttribute(*field, type, ghost_type);
        }
        *this << CloseTag{};
      }

      void dump(Support<Mesh> & support) override {
        auto && mesh_name = support.getName();
        auto && nodes = support.getNodes();

        if (support.hasProperty("xdmf_element") and
            support.getFields().empty()) {
          xi_include(support);
        } else {
          support.addProperty("xdmf_path", currentXMFPath());
          support.addProperty("xdmf_element", currentXMFElement());

          for (auto && type : support.elementTypes()) {
            /* clang-format off */
            *this << Tag("Grid")
                ("Name", mesh_name + ":" + std::to_string(type))
                ("GridType", "Uniform");
            /* clang-format on */

            // Topology
            auto & connectivities = support.getConnectivities().array(type);
            dumpTopology(connectivities, type);

            // Geometry
            dumpGeometry(nodes);
            for (auto && [_, field] : support.getFields()) {
              if (not field->is_visualization()) {
                continue;
              }

              if (field->is_nodal_field()) {
                this->dumpAttribute(*field);
              } else if (field->is_elemental_field()) {
                this->dumpAttribute(*field, type);
              }
            }

            for (auto && [_, sub_support] : support.getSubSupports()) {
              dumpSet(aka::as_type<Support<ElementGroup>>(*sub_support), type);
            }

            *this << CloseTag{};
          }
        }
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
