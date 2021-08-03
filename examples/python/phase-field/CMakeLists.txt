add_mesh(mesh_phasefield plate.geo DIM 2)
add_mesh(mesh_phasefield_static plate_static.geo DIM 2)

register_example(phasefield_dynamic
  SCRIPT phasefield-dynamic.py
  PYTHON
  FILES_TO_COPY material.dat
  DEPENDS mesh_phasefield
  )

register_example(phasefield_static
  SCRIPT phasefield-static.py
  PYTHON
  FILES_TO_COPY material_static.dat
  DEPENDS mesh_phasefield_static
  )
