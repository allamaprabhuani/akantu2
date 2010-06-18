/**
 * @file   vector.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jun 17 15:14:24 2010
 *
 * @brief  class vector
 *
 * @section LICENSE
 *
 * <insert lisence here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "common.hh"
#include "vector.hh"

/* -------------------------------------------------------------------------- */

void VectorBase::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(int i = 0; i < indent; i++, space += MYFEM_INDENT);
  stream << indent << "VectorBase [" << std::endl;
  stream << indent << " + nb tuples           : " << nb_tuples << std::endl;
  stream << indent << " + nb_component        : " << nb_component << std::endl;
  stream << indent << " + nb allocated tuples : " << nb_allocated_tuples << std::endl;
  int size = nb_allocated_tuples * nb_component * size_of_type;
  stream << indent << " + size : " << size << "B" << std::endl;
  stream << indent << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
