/**
 * @file   test_buildfacets_triangle_3.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 *
 * @date creation: Fri Sep 18 2015
 * @date last modification:  Wed Nov 08 2017
 *
 * @brief  Test to check the building of the facets. Mesh with triangles
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

template <typename types_>
class TestBuildFacetFixture : public ::testing::Test {
public:
  static constexpr const ElementType type1 =
      std::tuple_element_t<0, types_>::value;
  static constexpr const ElementType type2 =
      std::tuple_element_t<1, types_>::value;
  static constexpr const size_t dim =
      ElementClass<type1>::getSpatialDimension();

  void SetUp() override {
    mesh = std::make_unique<Mesh>(dim);
    mesh->read(std::to_string(type1) +
               ((type2 != _not_defined) ? std::to_string(type2) : "") + ".msh");
  }

  void TearDown() override { mesh.reset(nullptr); }

  template <class List>
  void print_connection(const std::string & name, ElementType type,
                        const List & list) {
    std::cout << name << std::endl;
    for (auto && outer : enumerate(list)) {
      std::cout << type << " " << std::get<0>(outer) << " connected to ";
      for (auto && inner : std::get<1>(outer)) {
        std::cout << inner << ", ";
      }
      std::cout << " " << std::endl;
    }
  }

protected:
  std::unique_ptr<Mesh> mesh;
};

template <typename type_>
constexpr const ElementType TestBuildFacetFixture<type_>::type1;
template <typename type_>
constexpr const ElementType TestBuildFacetFixture<type_>::type2;

template <typename type_>
constexpr const size_t TestBuildFacetFixture<type_>::dim;

template <typename T1, typename T2 = element_type_t<_not_defined>>
struct element_pair {
  using type = std::tuple<T1, T2>;
};

template <typename T1, typename T2 = element_type_t<_not_defined>>
using element_pair_t = typename element_pair<T1, T2>::type;

using buildfacets_types = ::testing::Types<
    element_pair_t<_element_type_triangle_3>,
    element_pair_t<_element_type_triangle_6>,
    element_pair_t<_element_type_quadrangle_4>,
    element_pair_t<_element_type_quadrangle_8>,
    element_pair_t<_element_type_tetrahedron_10>,
    element_pair_t<_element_type_hexahedron_8>,
    element_pair_t<_element_type_hexahedron_20>,
    element_pair_t<_element_type_pentahedron_6>,
    element_pair_t<_element_type_pentahedron_15>,
    element_pair_t<_element_type_quadrangle_4, _element_type_triangle_3>,
    element_pair_t<_element_type_quadrangle_8, _element_type_triangle_6>,
    element_pair_t<_element_type_hexahedron_8, _element_type_pentahedron_6>,
    element_pair_t<_element_type_hexahedron_20, _element_type_pentahedron_15>>;

TYPED_TEST_SUITE(TestBuildFacetFixture, buildfacets_types, );

TYPED_TEST(TestBuildFacetFixture, Buildfacets) {
  const auto & mesh_facets = this->mesh->initMeshFacets("mesh_facets");

  /* ------------------------------------------------------------------------ */
  /* Element to Subelement testing                                            */
  /* ------------------------------------------------------------------------ */
  const auto types_facet1{this->mesh->getAllFacetTypes(this->type1)};

  std::set<ElementType> types_facet;
  for (auto type : types_facet1) {
    types_facet.insert(type);
  }

  if (this->type2 != _not_defined) {
    const auto types_facet2{this->mesh->getAllFacetTypes(this->type2)};
    for (auto type : types_facet2) {
      types_facet.insert(type);
    }
  }

  for (auto type_facet : types_facet) {
    const auto & el_to_subel1 = mesh_facets.getElementToSubelement(type_facet);
    this->print_connection("ElementToSubelement1", type_facet, el_to_subel1);
  }

  const auto type_subfacet = this->mesh->getFacetType(*types_facet.begin());
  const auto & el_to_subel2 = mesh_facets.getElementToSubelement(type_subfacet);
  this->print_connection("ElementToSubelement2", type_subfacet, el_to_subel2);

  if (this->dim == 3) {
    const auto type_subsubfacet = this->mesh->getFacetType(type_subfacet);
    const auto & el_to_subel3 =
        mesh_facets.getElementToSubelement(type_subsubfacet);
    this->print_connection("ElementToSubelement3", type_subsubfacet,
                           el_to_subel3);
  }

  /* ------------------------------------------------------------------------ */
  /* Subelement to Element testing                                            */
  /* ------------------------------------------------------------------------ */
  std::cout << std::endl;
  const auto & subel_to_el1 = mesh_facets.getSubelementToElement(this->type1);
  this->print_connection("SubelementToElement1", this->type1,
                         MatrixProxy<Element>(subel_to_el1.data(),
                                              subel_to_el1.getNbComponent(),
                                              subel_to_el1.size()));

  if (this->type2 != _not_defined) {
    const auto & subel_to_el1 = mesh_facets.getSubelementToElement(this->type2);
    this->print_connection("SubelementToElement1", this->type2,
                           MatrixProxy<Element>(subel_to_el1.data(),
                                                subel_to_el1.getNbComponent(),
                                                subel_to_el1.size()));
  }

  for (auto type_facet : types_facet) {
    auto && subel_to_el2 = mesh_facets.getSubelementToElement(type_facet);
    this->print_connection("SubelementToElement2", type_facet,
                           MatrixProxy<Element>(subel_to_el2.data(),
                                                subel_to_el2.getNbComponent(),
                                                subel_to_el2.size()));
  }

  if (this->dim == 3) {
    const auto & subel_to_el3 =
        mesh_facets.getSubelementToElement(type_subfacet);
    this->print_connection("SubelementToElement3", type_subfacet,
                           MatrixProxy<Element>(subel_to_el3.data(),
                                                subel_to_el3.getNbComponent(),
                                                subel_to_el3.size()));
  }
}
