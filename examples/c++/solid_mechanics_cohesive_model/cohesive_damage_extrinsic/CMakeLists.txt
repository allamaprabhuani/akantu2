project(cohesive_damage_extrinsic)
cmake_minimum_required(VERSION 3.0.0)

find_package(Akantu REQUIRED)

add_akantu_simulation(cohesive_damage_extrinsic cohesive_damage_extrinsic.cc material.dat triangle.msh)
