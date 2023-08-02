/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
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
 */

/* -------------------------------------------------------------------------- */

#if defined(__INTEL_COMPILER)
//#pragma warning ( disable : 383 )
#elif defined(__clang__) // test clang to be sure that when we test for gnu it
// is only gnu
#elif (defined(__GNUC__) || defined(__GNUG__))
#define GCC_VERSION                                                            \
  (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION > 40600
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif

/* -------------------------------------------------------------------------- */
#include "parser.hh"
#if !defined(DOXYGEN)
#include "parser_grammar_tmpl.hh"
/* -------------------------------------------------------------------------- */
#include "algebraic_parser.hh"
#endif
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
RandomParameter<Real>
Parser::parseRandomParameter(const std::string & value,
                             const ParserSection & section) {
#if !defined(DOXYGEN)
  using boost::spirit::ascii::space_type;
  parser::RandomGeneratorGrammar<std::string::const_iterator, space_type>
      grammar(section);
  grammar.name("random_grammar");
  parser::ParsableRandomGenerator rg =
      Parser::parseType<parser::ParsableRandomGenerator>(value, grammar);

  Vector<Real> params = rg.parameters;
  switch (params.size()) {
  case 0:
    return make_random_parameter(rg.base, rg.type);
  case 1:
    return make_random_parameter(rg.base, rg.type, params(0));
  case 2:
    return make_random_parameter(rg.base, rg.type, params(0), params(1));
  }

  return make_random_parameter(rg.base, rg.type);
#endif
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
