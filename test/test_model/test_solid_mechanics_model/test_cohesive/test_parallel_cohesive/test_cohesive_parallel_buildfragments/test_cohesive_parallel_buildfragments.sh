#!/bin/bash

rm -r paraview
mkdir paraview

mpirun -np 4 ./test_cohesive_parallel_buildfragments
