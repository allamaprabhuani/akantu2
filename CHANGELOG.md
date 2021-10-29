## Version 4.0 (09-21-2021)

### Added
  
- pybind11 binding
- contact mechanics model
- phase field model
- Added a Changelog

### Changed

- transferred CI from jenkinsfile to gitlab CI/CD
- API changes to make container mode STL compatible
  - clear does not set to 0 anymore but empties containers
  - empty does not empty containers but tells if the container is empty
  - zero replace the old empty and set containers to 0
  
### Deprecated

- `getForce` in the `SolidMechanicsModel` becomes `getExternalForce`
- `firstType()`, `lastType()` replaced by `elementTypes()`

## Version 3.2 (not released)

### Added

- Activating PETSc solver back with the new solver interface

### Deprecated 

- deprecating old C++ 03 code


## 3.0 (2018-03)

### Added

- Parallel cohesive elements
- Element groups created by default for “physical_names”
- Named arguments for functions (e.g. model.initFull(_analysis_method = _static))

### Changed

- Models using new interface for solvers
  - Same configuration for all models
  - Solver can be configured in input file
- Only one function to solve a step model.solveStep()
- Simplification of the parallel simulation with the mesh.distribute() function
- Switch from C++ standard 2003 to 2014 Example of changes implied by this:

   for (UInt g = _not_ghost; g <= _ghost; ++g) {
      GhostType gt = (GhostType)g;
      Mesh::type_iterator it = this->mesh.firstType(spatial_dimension, gt);
      Mesh::type_iterator end = this->mesh.lastType(spatial_dimension, gt);
      for (; it != end; ++it) {
        ElementType & type = *it;
        ...
      }
    }

  becomes:

    for (auto ghost_type : ghost_types) {
      for (auto type : mesh.elementTypes(spatial_dimension,
                                         ghost_type)) {
        ...
      }
    }

### Deleted

- PETSc interface temporary inactive
- Periodic boundary condition temporary inactive

## 2.3 (2016-03)

### Added

- swig python interface

## 2.2 (2014-09)

### Added
- Cohesive elements

## 1.0 (2012-06)

### Added
- Continuum damage local and non-local
- Models: solid mechanics, structural mechanics, heat transfer

