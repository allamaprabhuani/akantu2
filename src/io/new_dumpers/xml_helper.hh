/**
 * @file   xml_helper.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Mon Oct 09 2017
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
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <tuple>
#include <type_traits>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_XML_HELPER_HH__
#define __AKANTU_XML_HELPER_HH__

namespace akantu {
namespace dumper {
  namespace XML {
    namespace {
      template <bool...> struct bool_pack;

      template <bool... v>
      using all_true =
          std::is_same<bool_pack<true, v...>, bool_pack<v..., true>>;

      template <std::size_t N> struct TupleUnpackHelper {
        template <typename... T>
        static void print(std::stringstream & sstr,
                          const std::tuple<T...> & t) {
          TupleUnpackHelper<N - 1>::print(sstr, t);
          sstr << " " << std::get<N - 1>(t);
        }
      };

      template <> struct TupleUnpackHelper<1> {
        template <typename... T>
        static void print(std::stringstream & sstr,
                          const std::tuple<T...> & t) {
          sstr << std::get<0>(t);
        }
      };

    } // namespace

    /* ---------------------------------------------------------------------- */
    class File;

    /* ---------------------------------------------------------------------- */
    class Tag {
    public:
      Tag(std::string tag, bool autoclosing_ = false)
          : tag_(std::move(tag)), autoclosing_(autoclosing_) {}

      Tag & operator()(const std::string & property,
                       const std::string_view & value) {
        properties.emplace(property, value);
        return *this;
      }

      template <typename T>
      Tag &
      operator()(const std::string & property, const T & value,
                 std::enable_if_t<
                     std::is_same<decltype(std::to_string(std::declval<T>())),
                                  std::string>::value> * /* unused */
                 = nullptr) {
        properties.emplace(property, std::to_string(value));
        return *this;
      }

      template <typename... T>
      Tag & operator()(
          const std::string & property, const std::tuple<T...> & value,
          std::enable_if_t<
              all_true<std::is_integral<T>::value...>::value> * /* unused */
          = nullptr) {
        std::stringstream sstr;
        TupleUnpackHelper<sizeof...(T)>::print(sstr, value);
        properties.emplace(property, sstr.str());
        return *this;
      }

      [[nodiscard]] const std::string &
      getProperty(const std::string & property) const {
        return properties.at(property);
      }

      [[nodiscard]] bool hasProperty(const std::string & property) const {
        return (properties.find(property) != properties.end());
      }

      [[nodiscard]] const std::string & tag() const { return tag_; }
      [[nodiscard]] bool autoclosing() const { return autoclosing_; }

    private:
      friend File;
      std::string tag_;
      bool autoclosing_{false};
      std::map<std::string, std::string> properties;
    };

    /* ---------------------------------------------------------------------- */
    struct CloseTag {};

    /* ---------------------------------------------------------------------- */
    class Indent {
    public:
      explicit Indent(char nb = 2) {
        for (char i = 0; i < nb; i++) {
          increment += " ";
        }
      }

      Indent & operator++() {
        indent += increment;
        return *this;
      }

      Indent & operator--() {
        indent = indent.substr(0, indent.size() - increment.size());
        return *this;
      }

      std::string str() const { return indent; }

    private:
      std::string increment{""};
      std::string indent{""};
    };

    /* ---------------------------------------------------------------------- */
    namespace fs = std::filesystem;

    class File {
    public:
      explicit File(const fs::path & path) {
        auto path_wof = path;
        path_wof.remove_filename();
        if (not fs::exists(path_wof)) {
          fs::create_directories(path_wof);
        }

        stream.open(path.string());
      }

      File(const File &) = delete;
      File(File &&) = delete;
      File & operator=(const File &) = delete;
      File & operator=(File &&) = delete;

      ~File() {
        auto size = tag_heap.size();
        for (std::size_t i = 0; i < size; ++i) {
          *this << CloseTag{};
        }

        stream.close();
      }

      virtual File & operator<<(const std::string & str) {
        if (new_line) {
          stream << indent.str();
          new_line = false;
        }

        auto pos_ret = str.find('\n');
        if (pos_ret == std::string::npos) { // not carry return
          stream << str;
        } else {
          stream << str.substr(0, pos_ret + 1);
          new_line = true;
          if (pos_ret != str.size() - 1) {
            this->operator<<(str.substr(pos_ret + 1));
          }
        }
        return *this;
      }

      virtual File & operator<<(const CloseTag & /* unused */) {
        --indent;
        *this << "</" << tag_heap.back() << ">\n";
        tag_heap.pop_back();
        return *this;
      }

      virtual File & operator<<(const Tag & t) {
        if (not t.autoclosing()) {
          tag_heap.push_back(t.tag());
        }
        *this << "<" << t.tag();

        std::for_each(t.properties.begin(), t.properties.end(), [&](auto && p) {
          *this << " " << p.first << "=\"" << p.second << "\"";
        });

        if (t.autoclosing()) {
          *this << " /";
        }

        *this << ">\n";

        if (not t.autoclosing()) {
          ++indent;
        }

        return *this;
      }

    private:
      std::ofstream stream;
      Indent indent{2};
      std::vector<std::string> tag_heap;
      bool new_line{true};
    };
  } // namespace XML
} // namespace dumper
} // namespace akantu

#endif /* __AKANTU_XML_HELPER_HH__ */

/* -------------------------------------------------------------------------- */
