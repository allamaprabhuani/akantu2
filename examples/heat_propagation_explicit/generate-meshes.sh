#! /bin/bash

for i in 16 32 64 128 192 256 320 396 456; do 
    export MESH_FILE=double_cube_tet4_$i.msh ; 
    export GEO_FILE=double_cube_tet_$i.geo ; 
    export GEO_FILE_TEMPLATE=double_cube_tet_template.geo ; 
    sed s/SSSSSS/$i/g $GEO_FILE_TEMPLATE > $GEO_FILE ;  
    gmsh -optimize -order 1 -3 $GEO_FILE -o $MESH_FILE ; 
done
