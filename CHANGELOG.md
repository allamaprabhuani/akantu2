## Version 5.0.3 (09-06-2023)

### Added

- Extending contact mechanics python interface

### Changed

- Bug fixes in conctact mechanics


## Version 5.0.2 (06-30-2023)

### Added

- pypi package for python 3.11 and 3.12

## Version 5.0.1 (05-05-2023)

### Changed

- Bug fixes:
  - bad constcorectness in arrays
  - applyBC broken

- Extending python API

## Version 5.0 (03-28-2023)

### Changed 

- C++ standard 17 is now accepted
- `Vector<T>` and `Matrix<T>` changed from internal types to Eigen::Matrix<T>
- `VectorProxy<T>` and `MatrixProxy<T>` are now `Eigen::Map<Eigen::Matrix<T>>`
  This introduces a potential bug if codes like the following example where used:
  ```
  for(auto && v_ : make_view(vectors, dim)) {
    Vector<Real> v(v_);
    ...
  }
  ```
  With the new version the temporary vector `v` will be a deep copy of `v_`
  instead of a shallow copy as in the previous version
  
### Added

- `make_view` as a static dimension version for vectors and matrices
  `make_view<size>(vectors)` and `make_view<rows, cols>(matrices)`
- `zip` iterators can be named, in which case the
  return tuple is a `named_tuple`.
  ```
  for(auto && t : zip("a"_n = as, "b"_n = bs)) {
    auto && a = t["a"_n];
    auto && b = t["b"_n];
    ...
  }
  ```
  
### Deprecated

- `begin_reinterpret` and `end_reinterpret` are error prone and deprecated in favor of `make_view`
- `storage()` members are deprecated in favor of `data()` in order to be compatible with the STL
- `get(.*)Energy(ElementType type, Idx index)` are deprecated in favor of `get(.*)Energy(const Element & element)`
  elements can be implicitly created, `getEnergy({type, index, _not_ghost})`
- `begin_(node|element)_group` and `end_(node|element)_group` are replaced by `iterate(Node|Element)Groups` 
- In the python interface:
  - the global `setDebugLevel`, `getDebugLevel` and `printBacktrace` were moved in the sub-module `debug`
  - the call to `finalize` is not needed
  - `applyDirichketBC` was replaced by `applyBC`

### Deleted

- `getForce`, `firstType()`, `lastType()` that where deprecated in version 4.0


## Version 4.0 (09-21-2021)

### Added
  
- pybind11 binding
- contact mechanics model
- phase field model
- Added a Changelog

### Changed

- transferred CI from jenkinsfile to gitlab CI/CD
- API changes to make container more STL compatible
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

   for (Int g = _not_ghost; g <= _ghost; ++g) {
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

