/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

// #include "element_class.hh"

// namespace akantu {

// /* --------------------------------------------------------------------------
// */
// template<> UInt ElementClass<_igfem_segment_3>::nb_sub_elements = 2;
// template<> UInt ElementClass<_igfem_triangle_4>::nb_sub_elements = 2;
// template<> UInt ElementClass<_igfem_triangle_5>::nb_sub_elements = 2;

// /* !!! stored as a matrix nb_subelements X nb_nodes_per_subelement in COL
// MAJOR */
// /* --------------------------------------------------------------------------
// */
// template<> UInt
// ElementClass<_igfem_segment_3>::sub_element_connectivity_vect[]  = { 0, //
// first type
// 										     2,
// 										     2, // second type
// 										     1};
// template<> UInt
// ElementClass<_igfem_triangle_4>::sub_element_connectivity_vect[]  = { 0, //
// first type
// 										     1,
// 										     3,
// 										     0, // second type
// 										     3,
// 										     2};
// template<> UInt
// ElementClass<_igfem_triangle_5>::sub_element_connectivity_vect[] = { 0, //
// first type
// 										    3,
// 										    4,
// 										    3, // second type
// 										    1,
// 										    2,
// 										    4};

// template<> UInt * ElementClass<_igfem_segment_3>::sub_element_connectivity[]
// = { &sub_element_connectivity_vect[0],
// 										 &sub_element_connectivity_vect[2] };
// template<> UInt * ElementClass<_igfem_triangle_4>::sub_element_connectivity[]
// = { &sub_element_connectivity_vect[0],
// 										 &sub_element_connectivity_vect[3] };
// template<> UInt * ElementClass<_igfem_triangle_5>::sub_element_connectivity[]
// = { &sub_element_connectivity_vect[0],
// 										 &sub_element_connectivity_vect[3] };
// /* --------------------------------------------------------------------------
// */

// } // namespace akantu
