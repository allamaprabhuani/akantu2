/* -------------------------------------------------------------------------- */
#include "vector_petsc.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_VECTOR_PETSC_TMPL_H_
#define __AKANTU_VECTOR_PETSC_TMPL_H_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class IDs, class Values>
void VectorPETSc::getValues(const IDs & idx, Values & values) const {
  if (idx.size() == 0)
    return;

  ISLocalToGlobalMapping is_ltog_map;
  PETSc_call(VecGetLocalToGlobalMapping, x, &is_ltog_map);

  PetscInt n;
  Array<PetscInt> lidx(idx.size());
  PETSc_call(ISGlobalToLocalMappingApply, is_ltog_map, IS_GTOLM_MASK,
             idx.size(), idx.storage(), &n, lidx.storage());

  getValuesLocal(lidx, values);
}

/* -------------------------------------------------------------------------- */
template <class IDs, class Values>
void VectorPETSc::getValuesLocal(const IDs & idx, Values & values) const {
  if (idx.size() == 0)
    return;

  Vec x_ghosted{nullptr};
  PETSc_call(VecGhostGetLocalForm, x, &x_ghosted);
  // VecScatterBegin(scatter, x, x_local, INSERT_VALUES, SCATTER_FORWARD);
  // VecScatterEnd(scatter, x, x_local, INSERT_VALUES, SCATTER_FORWARD);

  if (not x_ghosted) {
    const PetscScalar * array;
    PETSc_call(VecGetArrayRead, x, &array);

    for (auto && data : zip(idx, values)) {
      auto i = std::get<0>(data);
      if (i != -1) {
        std::get<1>(data) = array[i];
      }
    }

    PETSc_call(VecRestoreArrayRead, x, &array);
    return;
  }

  PETSc_call(VecSetOption, x_ghosted, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  PETSc_call(VecGetValues, x_ghosted, idx.size(), idx.storage(),
             values.storage());
  PETSc_call(VecGhostRestoreLocalForm, x, &x_ghosted);
}

/* -------------------------------------------------------------------------- */
template <class IDs, class Values>
void VectorPETSc::addValues(const IDs & gidx, const Values & values,
                            Real scale_factor) {
  Real * to_add = values.storage();
  Values scaled_array;
  if (scale_factor != 1.) {
    scaled_array = values;
    scaled_array *= scale_factor;
    to_add = scaled_array.storage();
  }

  PETSc_call(VecSetOption, x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  PETSc_call(VecSetValues, x, gidx.size(), gidx.storage(), to_add, ADD_VALUES);

  applyModifications();
}

/* -------------------------------------------------------------------------- */
template <class IDs, class Values>
void VectorPETSc::addValuesLocal(const IDs & lidx, const Values & values,
                                 Real scale_factor) {
  Vec x_ghosted{nullptr};
  PETSc_call(VecGhostGetLocalForm, x, &x_ghosted);

  if (not x_ghosted) {
    Real * to_add = values.storage();
    Array<Real> scaled_array;
    if (scale_factor != 1.) {
      scaled_array.copy(values, false);
      scaled_array *= scale_factor;
      to_add = scaled_array.storage();
    }

    PETSc_call(VecSetOption, x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
    PETSc_call(VecSetValuesLocal, x, lidx.size(), lidx.storage(), to_add,
               ADD_VALUES);
    return;
  }

  PETSc_call(VecGhostRestoreLocalForm, x, &x_ghosted);

  ISLocalToGlobalMapping is_ltog_map;
  PETSc_call(VecGetLocalToGlobalMapping, x, &is_ltog_map);

  Array<Int> gidx(lidx.size());
  PETSc_call(ISLocalToGlobalMappingApply, is_ltog_map, lidx.size(),
             lidx.storage(), gidx.storage());
  addValues(gidx, values, scale_factor);
}

} // namespace akantu

#endif // __AKANTU_VECTOR_PETSC_TMPL_H_
