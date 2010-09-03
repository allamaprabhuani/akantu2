/**
 * @file   test/static_memory.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jun 11 11:55:54 2010
 *
 * @brief  unit test for the StaticMemory class
 *
 */

/* -------------------------------------------------------------------------- */
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "aka_static_memory.hh"
#include "aka_vector.hh"

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  akantu::StaticMemory * st_mem = akantu::StaticMemory::getStaticMemory();

  akantu::Vector<int> & test_int = st_mem->smalloc<int>(0, "test_int", 1000, 3);

  test_int.resize(1050);

  test_int.resize(2000);

  std::cout << *st_mem << std::endl;

  st_mem->sfree(0, "test_int");

  exit(EXIT_SUCCESS);
}
