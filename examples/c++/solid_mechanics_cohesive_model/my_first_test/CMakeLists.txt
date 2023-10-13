project(my_first_test)
cmake_minimum_required(VERSION 3.0.0)

find_package(Akantu REQUIRED)

add_akantu_simulation(my_first_test my_first_test.cc material.dat triangle.msh)
