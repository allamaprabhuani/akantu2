add_mesh(compression_mesh ./mesh/compression.geo DIM 2 ORDER 1)

register_example(compression
  SCRIPT compression.py
  PYTHON
  FILES_TO_COPY ./materials/compression.dat
  DEPENDS compression_mesh)
