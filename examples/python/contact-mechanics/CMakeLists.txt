add_mesh(compression_mesh compression.geo DIM 2 ORDER 1)

register_example(compression
  SCRIPT compression.py
  PYTHON
  FILES_TO_COPY compression.dat
  DEPENDS compression_mesh)

register_example(detection_explicit.py
  SCRIPT detection-explicit.py
  PYTHON
  FILES_TO_COPY detection-explicit.dat
    detection-explicit.msh
  )
