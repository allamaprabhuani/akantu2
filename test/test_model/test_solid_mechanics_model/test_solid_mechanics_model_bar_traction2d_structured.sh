#!/bin/bash

rm energy.csv

./test_solid_mechanics_model_bar_traction2d_structured

if [ $? -eq 0 ]
then
    ./test_cst_energy.pl energy.csv 1e-3
else
    return $?
fi
