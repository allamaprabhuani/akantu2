/**
 * @file   mesh_io_msh.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Thu Jan 21 2016
 *
 * @brief  Read/Write for MSH files generated by gmsh
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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

/* -----------------------------------------------------------------------------
   Version (Legacy) 1.0

   $NOD
   number-of-nodes
   node-number x-coord y-coord z-coord
   ...
   $ENDNOD
   $ELM
   number-of-elements
   elm-number elm-type reg-phys reg-elem number-of-nodes node-number-list
   ...
   $ENDELM
   -----------------------------------------------------------------------------
   Version 2.1

   $MeshFormat
   version-number file-type data-size
   $EndMeshFormat
   $Nodes
   number-of-nodes
   node-number x-coord y-coord z-coord
   ...
   $EndNodes
   $Elements
   number-of-elements
   elm-number elm-type number-of-tags < tag > ... node-number-list
   ...
   $EndElements
   $PhysicalNames
   number-of-names
   physical-dimension physical-number "physical-name"
   ...
   $EndPhysicalNames
   $NodeData
   number-of-string-tags
   < "string-tag" >
   ...
   number-of-real-tags
   < real-tag >
   ...
   number-of-integer-tags
   < integer-tag >
   ...
   node-number value ...
   ...
   $EndNodeData
   $ElementData
   number-of-string-tags
   < "string-tag" >
   ...
   number-of-real-tags
   < real-tag >
   ...
   number-of-integer-tags
   < integer-tag >
   ...
   elm-number value ...
   ...
   $EndElementData
   $ElementNodeData
   number-of-string-tags
   < "string-tag" >
   ...
   number-of-real-tags
   < real-tag >
   ...
   number-of-integer-tags
   < integer-tag >
   ...
   elm-number number-of-nodes-per-element value ...
   ...
   $ElementEndNodeData

   -----------------------------------------------------------------------------
   elem-type

   1:  2-node line.
   2:  3-node triangle.
   3:  4-node quadrangle.
   4:  4-node tetrahedron.
   5:  8-node hexahedron.
   6:  6-node prism.
   7:  5-node pyramid.
   8:  3-node second order line
   9:  6-node second order triangle
   10: 9-node second order quadrangle
   11: 10-node second order tetrahedron
   12: 27-node second order hexahedron
   13: 18-node second order prism
   14: 14-node second order pyramid
   15: 1-node point.
   16: 8-node second order quadrangle
   17: 20-node second order hexahedron
   18: 15-node second order prism
   19: 13-node second order pyramid
   20: 9-node third order incomplete triangle
   21: 10-node third order triangle
   22: 12-node fourth order incomplete triangle
   23: 15-node fourth order triangle
   24: 15-node fifth order incomplete triangle
   25: 21-node fifth order complete triangle
   26: 4-node third order edge
   27: 5-node fourth order edge
   28: 6-node fifth order edge
   29: 20-node third order tetrahedron
   30: 35-node fourth order tetrahedron
   31: 56-node fifth order tetrahedron
   -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "mesh_io.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
// The boost spirit is a work on the way it does not compile so I kept the
// current code. The current code does not handle files generated on Windows
// <CRLF>

// #include <boost/config/warning_disable.hpp>
// #include <boost/spirit/include/qi.hpp>
// #include <boost/spirit/include/phoenix_core.hpp>
// #include <boost/spirit/include/phoenix_fusion.hpp>
// #include <boost/spirit/include/phoenix_object.hpp>
// #include <boost/spirit/include/phoenix_container.hpp>
// #include <boost/spirit/include/phoenix_operator.hpp>
// #include <boost/spirit/include/phoenix_bind.hpp>
// #include <boost/spirit/include/phoenix_stl.hpp>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
/*   Methods Implentations                                                    */
/* -------------------------------------------------------------------------- */
MeshIOMSH::MeshIOMSH() {
  canReadSurface = true;
  canReadExtendedData = true;

  _msh_nodes_per_elem[_msh_not_defined] = 0;
  _msh_nodes_per_elem[_msh_segment_2] = 2;
  _msh_nodes_per_elem[_msh_triangle_3] = 3;
  _msh_nodes_per_elem[_msh_quadrangle_4] = 4;
  _msh_nodes_per_elem[_msh_tetrahedron_4] = 4;
  _msh_nodes_per_elem[_msh_hexahedron_8] = 8;
  _msh_nodes_per_elem[_msh_prism_1] = 6;
  _msh_nodes_per_elem[_msh_pyramid_1] = 1;
  _msh_nodes_per_elem[_msh_segment_3] = 3;
  _msh_nodes_per_elem[_msh_triangle_6] = 6;
  _msh_nodes_per_elem[_msh_quadrangle_9] = 9;
  _msh_nodes_per_elem[_msh_tetrahedron_10] = 10;
  _msh_nodes_per_elem[_msh_hexahedron_27] = 27;
  _msh_nodes_per_elem[_msh_hexahedron_20] = 20;
  _msh_nodes_per_elem[_msh_prism_18] = 18;
  _msh_nodes_per_elem[_msh_prism_15] = 15;
  _msh_nodes_per_elem[_msh_pyramid_14] = 14;
  _msh_nodes_per_elem[_msh_point] = 1;
  _msh_nodes_per_elem[_msh_quadrangle_8] = 8;

  _msh_to_akantu_element_types[_msh_not_defined] = _not_defined;
  _msh_to_akantu_element_types[_msh_segment_2] = _segment_2;
  _msh_to_akantu_element_types[_msh_triangle_3] = _triangle_3;
  _msh_to_akantu_element_types[_msh_quadrangle_4] = _quadrangle_4;
  _msh_to_akantu_element_types[_msh_tetrahedron_4] = _tetrahedron_4;
  _msh_to_akantu_element_types[_msh_hexahedron_8] = _hexahedron_8;
  _msh_to_akantu_element_types[_msh_prism_1] = _pentahedron_6;
  _msh_to_akantu_element_types[_msh_pyramid_1] = _not_defined;
  _msh_to_akantu_element_types[_msh_segment_3] = _segment_3;
  _msh_to_akantu_element_types[_msh_triangle_6] = _triangle_6;
  _msh_to_akantu_element_types[_msh_quadrangle_9] = _not_defined;
  _msh_to_akantu_element_types[_msh_tetrahedron_10] = _tetrahedron_10;
  _msh_to_akantu_element_types[_msh_hexahedron_27] = _not_defined;
  _msh_to_akantu_element_types[_msh_hexahedron_20] = _hexahedron_20;
  _msh_to_akantu_element_types[_msh_prism_18] = _not_defined;
  _msh_to_akantu_element_types[_msh_prism_15] = _pentahedron_15;
  _msh_to_akantu_element_types[_msh_pyramid_14] = _not_defined;
  _msh_to_akantu_element_types[_msh_point] = _point_1;
  _msh_to_akantu_element_types[_msh_quadrangle_8] = _quadrangle_8;

  _akantu_to_msh_element_types[_not_defined] = _msh_not_defined;
  _akantu_to_msh_element_types[_segment_2] = _msh_segment_2;
  _akantu_to_msh_element_types[_segment_3] = _msh_segment_3;
  _akantu_to_msh_element_types[_triangle_3] = _msh_triangle_3;
  _akantu_to_msh_element_types[_triangle_6] = _msh_triangle_6;
  _akantu_to_msh_element_types[_tetrahedron_4] = _msh_tetrahedron_4;
  _akantu_to_msh_element_types[_tetrahedron_10] = _msh_tetrahedron_10;
  _akantu_to_msh_element_types[_quadrangle_4] = _msh_quadrangle_4;
  _akantu_to_msh_element_types[_quadrangle_8] = _msh_quadrangle_8;
  _akantu_to_msh_element_types[_hexahedron_8] = _msh_hexahedron_8;
  _akantu_to_msh_element_types[_hexahedron_20] = _msh_hexahedron_20;
  _akantu_to_msh_element_types[_pentahedron_6] = _msh_prism_1;
  _akantu_to_msh_element_types[_pentahedron_15] = _msh_prism_15;
  _akantu_to_msh_element_types[_point_1] = _msh_point;

#if defined(AKANTU_STRUCTURAL_MECHANICS)
  _akantu_to_msh_element_types[_bernoulli_beam_2] = _msh_segment_2;
  _akantu_to_msh_element_types[_bernoulli_beam_3] = _msh_segment_2;
  _akantu_to_msh_element_types[_kirchhoff_shell] = _msh_triangle_3;
#endif

  std::map<ElementType, MSHElementType>::iterator it;
  for (it = _akantu_to_msh_element_types.begin();
       it != _akantu_to_msh_element_types.end(); ++it) {
    UInt nb_nodes = _msh_nodes_per_elem[it->second];

    std::vector<UInt> tmp(nb_nodes);
    for (UInt i = 0; i < nb_nodes; ++i) {
      tmp[i] = i;
    }

    switch (it->first) {
    case _tetrahedron_10:
      tmp[8] = 9;
      tmp[9] = 8;
      break;
    case _pentahedron_6:
      tmp[0] = 2;
      tmp[1] = 0;
      tmp[2] = 1;
      tmp[3] = 5;
      tmp[4] = 3;
      tmp[5] = 4;
      break;
    case _pentahedron_15:
      tmp[0] = 2;
      tmp[1] = 0;
      tmp[2] = 1;
      tmp[3] = 5;
      tmp[4] = 3;
      tmp[5] = 4;
      tmp[6] = 8;
      tmp[8] = 11;
      tmp[9] = 6;
      tmp[10] = 9;
      tmp[11] = 10;
      tmp[12] = 14;
      tmp[14] = 12;
      break;
    case _hexahedron_20:
      tmp[9] = 11;
      tmp[10] = 12;
      tmp[11] = 9;
      tmp[12] = 13;
      tmp[13] = 10;
      tmp[17] = 19;
      tmp[18] = 17;
      tmp[19] = 18;
      break;
    default:
      // nothing to change
      break;
    }
    _read_order[it->first] = tmp;
  }
}

/* -------------------------------------------------------------------------- */
MeshIOMSH::~MeshIOMSH() {}

/* -------------------------------------------------------------------------- */
/* Spirit stuff                                                               */
/* -------------------------------------------------------------------------- */
// namespace _parser {
//   namespace spirit  = ::boost::spirit;
//   namespace qi      = ::boost::spirit::qi;
//   namespace ascii   = ::boost::spirit::ascii;
//   namespace lbs     = ::boost::spirit::qi::labels;
//   namespace phx     = ::boost::phoenix;

//   /* ------------------------------------------------------------------------
//   */
//   /* Lazy functors */
//   /* ------------------------------------------------------------------------
//   */
//   struct _Element {
//     int index;
//     std::vector<int> tags;
//     std::vector<int> connectivity;
//     ElementType type;
//   };

//   /* ------------------------------------------------------------------------
//   */
//   struct lazy_get_nb_nodes_ {
//     template <class elem_type> struct result { typedef int type; };
//     template <class elem_type> bool operator()(elem_type et) {
//       return MeshIOMSH::_msh_nodes_per_elem[et];
//     }
//   };

//   /* ------------------------------------------------------------------------
//   */
//   struct lazy_element_ {
//     template <class id_t, class tags_t, class elem_type, class conn_t>
//     struct result {
//       typedef _Element type;
//     };
//     template <class id_t, class tags_t, class elem_type, class conn_t>
//     _Element operator()(id_t id, const elem_type & et, const tags_t & t,
//                        const conn_t & c) {
//       _Element tmp_el;
//       tmp_el.index = id;
//       tmp_el.tags = t;
//       tmp_el.connectivity = c;
//       tmp_el.type = et;
//       return tmp_el;
//     }
//   };

//   /* ------------------------------------------------------------------------
//   */
//   struct lazy_check_value_ {
//     template <class T> struct result { typedef void type; };
//     template <class T> void operator()(T v1, T v2) {
//       if (v1 != v2) {
//         AKANTU_EXCEPTION("The msh parser expected a "
//                          << v2 << " in the header bug got a " << v1);
//       }
//     }
//   };

//   /* ------------------------------------------------------------------------
//   */
//   struct lazy_node_read_ {
//     template <class Mesh, class ID, class V, class size, class Map>
//     struct result {
//       typedef bool type;
//     };
//     template <class Mesh, class ID, class V, class size, class Map>
//     bool operator()(Mesh & mesh, const ID & id, const V & pos, size max,
//                     Map & nodes_mapping) const {
//       Vector<Real> tmp_pos(mesh.getSpatialDimension());
//       UInt i = 0;
//       for (typename V::const_iterator it = pos.begin();
//            it != pos.end() || i < mesh.getSpatialDimension(); ++it)
//         tmp_pos[i++] = *it;

//       nodes_mapping[id] = mesh.getNbNodes();
//       mesh.getNodes().push_back(tmp_pos);
//       return (mesh.getNbNodes() < UInt(max));
//     }
//   };

//   /* ------------------------------------------------------------------------
//   */
//   struct lazy_element_read_ {
//     template <class Mesh, class EL, class size, class NodeMap, class ElemMap>
//     struct result {
//       typedef bool type;
//     };

//     template <class Mesh, class EL, class size, class NodeMap, class ElemMap>
//     bool operator()(Mesh & mesh, const EL & element, size max,
//                     const NodeMap & nodes_mapping,
//                     ElemMap & elements_mapping) const {
//       Vector<UInt> tmp_conn(Mesh::getNbNodesPerElement(element.type));

//       AKANTU_DEBUG_ASSERT(element.connectivity.size() == tmp_conn.size(),
//                           "The element "
//                               << element.index
//                               << "in the MSH file has too many nodes.");

//       mesh.addConnectivityType(element.type);
//       Array<UInt> & connectivity = mesh.getConnectivity(element.type);

//       UInt i = 0;
//       for (std::vector<int>::const_iterator it =
//       element.connectivity.begin();
//            it != element.connectivity.end(); ++it) {
//         typename NodeMap::const_iterator nit = nodes_mapping.find(*it);
//         AKANTU_DEBUG_ASSERT(nit != nodes_mapping.end(),
//                             "There is an unknown node in the connectivity.");
//         tmp_conn[i++] = nit->second;
//       }

//       ::akantu::Element el(element.type, connectivity.size());
//       elements_mapping[element.index] = el;

//       connectivity.push_back(tmp_conn);

//       for (UInt i = 0; i < element.tags.size(); ++i) {
//         std::stringstream tag_sstr;
//         tag_sstr << "tag_" << i;
//         Array<UInt> * data =
//             mesh.template getDataPointer<UInt>(tag_sstr.str(), element.type,
//             _not_ghost);
//         data->push_back(element.tags[i]);
//       }

//       return (mesh.getNbElement() < UInt(max));
//     }
//   };

//   /* ------------------------------------------------------------------------
//   */
//   template <class Iterator, typename Skipper = ascii::space_type>
//   struct MshMeshGrammar : qi::grammar<Iterator, void(), Skipper> {
//   public:
//     MshMeshGrammar(Mesh & mesh)
//         : MshMeshGrammar::base_type(start, "msh_mesh_reader"), mesh(mesh) {
//       phx::function<lazy_element_> lazy_element;
//       phx::function<lazy_get_nb_nodes_> lazy_get_nb_nodes;
//       phx::function<lazy_check_value_> lazy_check_value;
//       phx::function<lazy_node_read_> lazy_node_read;
//       phx::function<lazy_element_read_> lazy_element_read;

// clang-format off

//       start
//         =  *( known_section | unknown_section
//             )
//        ;

//       known_section
//         =  qi::char_("$")
//            >> sections            [ lbs::_a = lbs::_1 ]
//            >> qi::lazy(*lbs::_a)
//            >> qi::lit("$End")
//           //>> qi::lit(phx::ref(lbs::_a))
//         ;

//       unknown_section
//         =  qi::char_("$")
//            >> qi::char_("")     [ lbs::_a = lbs::_1 ]
//            >> ignore_section
//            >> qi::lit("$End")
//            >> qi::lit(phx::val(lbs::_a))
//         ;

//       mesh_format // version followed by 0 or 1 for ascii or binary
//         =  version >> (
//                         ( qi::char_("0")
//                           >> qi::int_       [ lazy_check_value(lbs::_1, 8) ]
//                         )
//                         | ( qi::char_("1")
//                             >> qi::int_     [ lazy_check_value(lbs::_1, 8) ]
//                             >> qi::dword    [ lazy_check_value(lbs::_1, 1) ]
//                           )
//                         )
//         ;

//       nodes
//         =  nodes_
//         ;

//       nodes_
//         =  qi::int_                      [ lbs::_a = lbs::_1 ]
//            > *(
//                 ( qi::int_ >> position ) [ lazy_node_read(phx::ref(mesh),
//                                                           lbs::_1,
//                                                           phx::cref(lbs::_2),
//                                                           lbs::_a,
//                                                           phx::ref(this->msh_nodes_to_akantu)) ]
//               )
//         ;

//       element
//         =  elements_
//         ;

//       elements_
//         =  qi::int_                                 [ lbs::_a = lbs::_1 ]
//            > qi::repeat(phx::ref(lbs::_a))[ element [ lazy_element_read(phx::ref(mesh),
//                                                                         lbs::_1,
//                                                                         phx::cref(lbs::_a),
//                                                                         phx::cref(this->msh_nodes_to_akantu),
//                                                                         phx::ref(this->msh_elemt_to_akantu)) ]]
//         ;

//       ignore_section
//         =  *(qi::char_ - qi::char_('$'))
//         ;

//       interpolation_scheme = ignore_section;
//       element_data         = ignore_section;
//       node_data            = ignore_section;

//       version
//         =  qi::int_                         [ phx::push_back(lbs::_val, lbs::_1) ]
//            >> *( qi::char_(".") >> qi::int_ [ phx::push_back(lbs::_val, lbs::_1) ] )
//         ;

//       position
//         =  real   [ phx::push_back(lbs::_val, lbs::_1) ]
//            > real [ phx::push_back(lbs::_val, lbs::_1) ]
//            > real [ phx::push_back(lbs::_val, lbs::_1) ]
//         ;

//       tags
//         =  qi::int_      [ lbs::_a = lbs::_1 ]
//           > qi::repeat(phx::val(lbs::_a))[ qi::int_ [ phx::push_back(lbs::_val,
//                                                                      lbs::_1) ] ]
//         ;

//       element
//         =  ( qi::int_                [ lbs::_a = lbs::_1 ]
//              > msh_element_type      [ lbs::_b = lazy_get_nb_nodes(lbs::_1) ]
//              > tags                  [ lbs::_c = lbs::_1 ]
//              > connectivity(lbs::_a) [ lbs::_d = lbs::_1 ]
//            )                         [ lbs::_val = lazy_element(lbs::_a,
//                                                                 phx::cref(lbs::_b),
//                                                                 phx::cref(lbs::_c),
//                                                                 phx::cref(lbs::_d)) ]
//         ;
//       connectivity
//         =  qi::repeat(lbs::_r1)[ qi::int_ [ phx::push_back(lbs::_val,
//                                                            lbs::_1) ] ]
//         ;

//       sections.add
//           ("MeshFormat", &mesh_format)
//           ("Nodes", &nodes)
//           ("Elements", &elements)
//           ("PhysicalNames", &physical_names)
//           ("InterpolationScheme", &interpolation_scheme)
//           ("ElementData", &element_data)
//           ("NodeData", &node_data);

//       msh_element_type.add
//           ("0" , _not_defined   )
//           ("1" , _segment_2     )
//           ("2" , _triangle_3    )
//           ("3" , _quadrangle_4  )
//           ("4" , _tetrahedron_4 )
//           ("5" , _hexahedron_8  )
//           ("6" , _pentahedron_6 )
//           ("7" , _not_defined   )
//           ("8" , _segment_3     )
//           ("9" , _triangle_6    )
//           ("10", _not_defined   )
//           ("11", _tetrahedron_10)
//           ("12", _not_defined   )
//           ("13", _not_defined   )
//           ("14", _hexahedron_20 )
//           ("15", _pentahedron_15)
//           ("16", _not_defined   )
//           ("17", _point_1       )
//           ("18", _quadrangle_8  );

//       mesh_format         .name("MeshFormat"         );
//       nodes               .name("Nodes"              );
//       elements            .name("Elements"           );
//       physical_names      .name("PhysicalNames"      );
//       interpolation_scheme.name("InterpolationScheme");
//       element_data        .name("ElementData"        );
//       node_data           .name("NodeData"           );

// clang-format on
//     }

//     /* ----------------------------------------------------------------------
//     */
//     /* Rules */
//     /* ----------------------------------------------------------------------
//     */
//   private:
//     qi::symbols<char, ElementType> msh_element_type;
//     qi::symbols<char, qi::rule<Iterator, void(), Skipper> *> sections;
//     qi::rule<Iterator, void(), Skipper> start;
//     qi::rule<Iterator, void(), Skipper, qi::locals<std::string> >
//     unknown_section;
//     qi::rule<Iterator, void(), Skipper, qi::locals<qi::rule<Iterator,
//     Skipper> *> > known_section;
//     qi::rule<Iterator, void(), Skipper> mesh_format, nodes, elements,
//     physical_names, ignore_section,
//         interpolation_scheme, element_data, node_data, any_line;
//     qi::rule<Iterator, void(), Skipper, qi::locals<int> > nodes_;
//     qi::rule<Iterator, void(), Skipper, qi::locals< int, int, vector<int>,
//     vector<int> > > elements_;

//     qi::rule<Iterator, std::vector<int>(), Skipper> version;
//     qi::rule<Iterator, _Element(), Skipper, qi::locals<ElementType> >
//     element;
//     qi::rule<Iterator, std::vector<int>(), Skipper, qi::locals<int> > tags;
//     qi::rule<Iterator, std::vector<int>(int), Skipper> connectivity;
//     qi::rule<Iterator, std::vector<Real>(), Skipper> position;

//     qi::real_parser<Real, qi::real_policies<Real> > real;

//     /* ----------------------------------------------------------------------
//     */
//     /* Members */
//     /* ----------------------------------------------------------------------
//     */
//   private:
//     /// reference to the mesh to read
//     Mesh & mesh;

//     /// correspondance between the numbering of nodes in the abaqus file and
//     in
//     /// the akantu mesh
//     std::map<UInt, UInt> msh_nodes_to_akantu;

//     /// correspondance between the element number in the abaqus file and the
//     /// Element in the akantu mesh
//     std::map<UInt, Element> msh_elemt_to_akantu;
//   };
// }

// /* --------------------------------------------------------------------------
// */
// void MeshIOAbaqus::read(const std::string& filename, Mesh& mesh) {
//   namespace qi      = boost::spirit::qi;
//   namespace ascii   = boost::spirit::ascii;

//   std::ifstream infile;
//   infile.open(filename.c_str());

//   if(!infile.good()) {
//     AKANTU_DEBUG_ERROR("Cannot open file " << filename);
//   }

//   std::string storage; // We will read the contents here.
//   infile.unsetf(std::ios::skipws); // No white space skipping!
//   std::copy(std::istream_iterator<char>(infile),
//             std::istream_iterator<char>(),
//             std::back_inserter(storage));

//   typedef std::string::const_iterator iterator_t;
//   typedef ascii::space_type skipper;
//   typedef _parser::MshMeshGrammar<iterator_t, skipper> grammar;

//   grammar g(mesh);
//   skipper ws;

//   iterator_t iter = storage.begin();
//   iterator_t end  = storage.end();

//   qi::phrase_parse(iter, end, g, ws);

//   this->setNbGlobalNodes(mesh, mesh.getNodes().size());
//   MeshUtils::fillElementToSubElementsData(mesh);
// }

static void my_getline(std::ifstream & infile, std::string & str) {
  std::string tmp_str;
  std::getline(infile, tmp_str);
  str = trim(tmp_str);
}

/* -------------------------------------------------------------------------- */
void MeshIOMSH::read(const std::string & filename, Mesh & mesh) {

  MeshAccessor mesh_accessor(mesh);

  std::ifstream infile;
  infile.open(filename.c_str());

  std::string line;
  UInt first_node_number = std::numeric_limits<UInt>::max(),
       last_node_number = 0, file_format = 1, current_line = 0;
  bool has_physical_names = false;

  if (!infile.good()) {
    AKANTU_DEBUG_ERROR("Cannot open file " << filename);
  }

  while (infile.good()) {
    my_getline(infile, line);
    current_line++;

    /// read the header
    if (line == "$MeshFormat") {
      my_getline(infile, line); /// the format line
      std::stringstream sstr(line);
      std::string version;
      sstr >> version;
      Int format;
      sstr >> format;
      if (format != 0)
        AKANTU_DEBUG_ERROR("This reader can only read ASCII files.");
      my_getline(infile, line); /// the end of block line
      current_line += 2;
      file_format = 2;
    }

    /// read the physical names
    if (line == "$PhysicalNames") {
      has_physical_names = true;
      my_getline(infile, line); /// the format line
      std::stringstream sstr(line);

      UInt num_of_phys_names;
      sstr >> num_of_phys_names;

      for (UInt k(0); k < num_of_phys_names; k++) {
        my_getline(infile, line);
        std::stringstream sstr_phys_name(line);
        UInt phys_name_id;
        UInt phys_dim;

        sstr_phys_name >> phys_dim >> phys_name_id;

        std::size_t b = line.find("\"");
        std::size_t e = line.rfind("\"");
        std::string phys_name = line.substr(b + 1, e - b - 1);

        phys_name_map[phys_name_id] = phys_name;
      }
    }

    /// read all nodes
    if (line == "$Nodes" || line == "$NOD") {
      UInt nb_nodes;

      my_getline(infile, line);
      std::stringstream sstr(line);
      sstr >> nb_nodes;
      current_line++;

      Array<Real> & nodes = mesh_accessor.getNodes();
      nodes.resize(nb_nodes);
      mesh_accessor.setNbGlobalNodes(nb_nodes);

      UInt index;
      Real coord[3];
      UInt spatial_dimension = nodes.getNbComponent();
      /// for each node, read the coordinates
      for (UInt i = 0; i < nb_nodes; ++i) {
        UInt offset = i * spatial_dimension;

        my_getline(infile, line);
        std::stringstream sstr_node(line);
        sstr_node >> index >> coord[0] >> coord[1] >> coord[2];
        current_line++;

        first_node_number = std::min(first_node_number, index);
        last_node_number = std::max(last_node_number, index);

        /// read the coordinates
        for (UInt j = 0; j < spatial_dimension; ++j)
          nodes.storage()[offset + j] = coord[j];
      }
      my_getline(infile, line); /// the end of block line
    }

    /// read all elements
    if (line == "$Elements" || line == "$ELM") {
      UInt nb_elements;

      std::vector<UInt> read_order;

      my_getline(infile, line);
      std::stringstream sstr(line);
      sstr >> nb_elements;
      current_line++;

      Int index;
      UInt msh_type;
      ElementType akantu_type, akantu_type_old = _not_defined;
      Array<UInt> * connectivity = NULL;
      UInt node_per_element = 0;

      for (UInt i = 0; i < nb_elements; ++i) {
        my_getline(infile, line);
        std::stringstream sstr_elem(line);
        current_line++;

        sstr_elem >> index;
        sstr_elem >> msh_type;

        /// get the connectivity vector depending on the element type
        akantu_type = this->_msh_to_akantu_element_types[(MSHElementType)msh_type];

        if (akantu_type == _not_defined) {
          AKANTU_DEBUG_WARNING("Unsuported element kind "
                               << msh_type << " at line " << current_line);
          continue;
        }

        if (akantu_type != akantu_type_old) {
          connectivity = &mesh_accessor.getConnectivity(akantu_type);
          //          connectivity->resize(0);

          node_per_element = connectivity->getNbComponent();
          akantu_type_old = akantu_type;
          read_order = this->_read_order[akantu_type];
        }

        /// read tags informations
        if (file_format == 2) {
          UInt nb_tags;
          sstr_elem >> nb_tags;
          for (UInt j = 0; j < nb_tags; ++j) {
            Int tag;
            sstr_elem >> tag;
            std::stringstream sstr_tag_name;
            sstr_tag_name << "tag_" << j;
            Array<UInt> & data = mesh.getDataPointer<UInt>(
                sstr_tag_name.str(), akantu_type, _not_ghost);
            data.push_back(tag);
          }
        } else if (file_format == 1) {
          Int tag;
          sstr_elem >> tag; // reg-phys
          std::string tag_name = "tag_0";
          Array<UInt> * data =
              &mesh.getDataPointer<UInt>(tag_name, akantu_type, _not_ghost);
          data->push_back(tag);

          sstr_elem >> tag; // reg-elem
          tag_name = "tag_1";
          data = &mesh.getDataPointer<UInt>(tag_name, akantu_type, _not_ghost);
          data->push_back(tag);

          sstr_elem >> tag; // number-of-nodes
        }

        Vector<UInt> local_connect(node_per_element);
        for (UInt j = 0; j < node_per_element; ++j) {
          UInt node_index;
          sstr_elem >> node_index;

          AKANTU_DEBUG_ASSERT(node_index <= last_node_number,
                              "Node number not in range : line "
                              << current_line);

          node_index -= first_node_number;
          local_connect(read_order[j]) = node_index;
        }
        connectivity->push_back(local_connect);
      }
      my_getline(infile, line); /// the end of block line
    }

    if ((line[0] == '$') && (line.find("End") == std::string::npos)) {
      AKANTU_DEBUG_WARNING("Unsuported block_kind " << line << " at line "
                                                    << current_line);
    }
  }

  //mesh.updateTypesOffsets(_not_ghost);

  infile.close();

  this->constructPhysicalNames("tag_0", mesh);

  if (has_physical_names)
    mesh.createGroupsFromMeshData<std::string>("physical_names");

  MeshUtils::fillElementToSubElementsData(mesh);
}

/* -------------------------------------------------------------------------- */
void MeshIOMSH::write(const std::string & filename, const Mesh & mesh) {
  std::ofstream outfile;
  const Array<Real> & nodes = mesh.getNodes();

  outfile.open(filename.c_str());

  outfile << "$MeshFormat" << std::endl;
  outfile << "2.1 0 8" << std::endl;
  outfile << "$EndMeshFormat" << std::endl;

  outfile << std::setprecision(std::numeric_limits<Real>::digits10);
  outfile << "$Nodes" << std::endl;

  outfile << nodes.size() << std::endl;

  outfile << std::uppercase;
  for (UInt i = 0; i < nodes.size(); ++i) {
    Int offset = i * nodes.getNbComponent();
    outfile << i + 1;
    for (UInt j = 0; j < nodes.getNbComponent(); ++j) {
      outfile << " " << nodes.storage()[offset + j];
    }

    for (UInt p = nodes.getNbComponent(); p < 3; ++p)
      outfile << " " << 0.;
    outfile << std::endl;
    ;
  }
  outfile << std::nouppercase;
  outfile << "$EndNodes" << std::endl;
  ;

  outfile << "$Elements" << std::endl;
  ;

  Mesh::type_iterator it =
      mesh.firstType(_all_dimensions, _not_ghost, _ek_not_defined);
  Mesh::type_iterator end =
      mesh.lastType(_all_dimensions, _not_ghost, _ek_not_defined);

  Int nb_elements = 0;
  for (; it != end; ++it) {
    const Array<UInt> & connectivity = mesh.getConnectivity(*it, _not_ghost);
    nb_elements += connectivity.size();
  }
  outfile << nb_elements << std::endl;

  UInt element_idx = 1;
  for (it = mesh.firstType(_all_dimensions, _not_ghost, _ek_not_defined);
       it != end; ++it) {
    ElementType type = *it;
    const Array<UInt> & connectivity = mesh.getConnectivity(type, _not_ghost);

    UInt * tag[2] = {NULL, NULL};
    try {
      const Array<UInt> & data_tag_0 =
          mesh.getData<UInt>("tag_0", type, _not_ghost);
      tag[0] = data_tag_0.storage();
    } catch (...) {
      tag[0] = NULL;
    }

    try {
      const Array<UInt> & data_tag_1 =
          mesh.getData<UInt>("tag_1", type, _not_ghost);
      tag[1] = data_tag_1.storage();
    } catch (...) {
      tag[1] = NULL;
    }

    for (UInt i = 0; i < connectivity.size(); ++i) {
      UInt offset = i * connectivity.getNbComponent();
      outfile << element_idx << " " << _akantu_to_msh_element_types[type]
              << " 2";

      /// \todo write the real data in the file
      for (UInt t = 0; t < 2; ++t)
        if (tag[t])
          outfile << " " << tag[t][i];
        else
          outfile << " 0";

      for (UInt j = 0; j < connectivity.getNbComponent(); ++j) {
        outfile << " " << connectivity.storage()[offset + j] + 1;
      }
      outfile << std::endl;
      element_idx++;
    }
  }

  outfile << "$EndElements" << std::endl;
  ;

  outfile.close();
}

/* -------------------------------------------------------------------------- */

} // akantu
