/**
 * @file   test_static_memory.cpp
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jun 11 11:55:54 2010
 *
 * @brief  unit test for the StaticMemory class
 *
 */

/* -------------------------------------------------------------------------- */
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "static_memory.hpp"

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  myfem::StaticMemory * st_mem = myfem::StaticMemory::getStaticMemory();

  int * test_int = st_mem->smalloc_int(0, "test_int", 1000, 3);



}
