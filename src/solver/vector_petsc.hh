/**
 * @file   vector_petsc.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Jan 01 2019
 *
 * @brief A Documented file.
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include "mpi_communicator_data.hh"
#include "petsc_wrapper.hh"
/* -------------------------------------------------------------------------- */
#include <petscvec.h>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_VECTOR_PETSC_HH__
#define __AKANTU_VECTOR_PETSC_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
namespace internal {
  /* ------------------------------------------------------------------------ */
  class VectorPETSc {
  public:
    virtual ~VectorPETSc() = default;

    operator Vec &() { return x; }
    operator const Vec &() const { return x; }

    virtual Int size() const {
      PetscInt n;
      PETSc_call(VecGetSize, x, &n);
      return n;
    }

    virtual Int local_size() const {
      PetscInt n;
      PETSc_call(VecGetLocalSize, x, &n);
      return n;
    }

    void getLocalToGlobalMapping(ISLocalToGlobalMapping & is_ltog_map) {
      PETSc_call(VecGetLocalToGlobalMapping, x, &is_ltog_map);
    }

    AKANTU_GET_MACRO_NOT_CONST(Vec, x, auto &);
    AKANTU_GET_MACRO(Vec, x, const auto &);

    virtual void printself(std::ostream & stream, int indent = 0) const {
      this->print(x, stream, indent);
    }

  protected:
    void print(const Vec & y, std::ostream & stream, int indent = 0) const;

  protected:
    Vec x{nullptr};
  };

} // namespace internal

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
class VectorPETSc : public internal::VectorPETSc {
public:
  VectorPETSc(const Communicator & communicator, const UInt n = 0,
              const SizeType & size_type = SizeType::_local,
              const ID & id = "vector_petsc");

  VectorPETSc(const VectorPETSc & Vector,
              const ID & id = "solver_vector_petsc");

  ~VectorPETSc() override;

  // resize the Vector to the size of the problem
  // void resize() override;
  void clear();

  VectorPETSc & operator=(const VectorPETSc & y);
  VectorPETSc & operator+(const VectorPETSc & y);

  /// get values using processors global indexes
  template <class IDs, class Values>
  void getValues(const IDs & idx, Values & values) const;

  /// get values using processors local indexes
  template <class IDs, class Values>
  void getValuesLocal(const IDs & idx, Values & values) const;

  /// adding values to the Vector using the global indices
  template <class IDs, class Values>
  void addValues(const IDs & gidx, const Values & values,
                 Real scale_factor = 1.);

  /// adding values to the Vector using the local indices
  template <class IDs, class Values>
  void addValuesLocal(const IDs & lidx, const Values & values,
                      Real scale_factor = 1.);

  /* ------------------------------------------------------------------------ */
  operator const Array<Real> &() const;

protected:
  void applyModifications();
  void updateGhost();

protected:
  MPI_Comm mpi_comm;
  Int n{0}, n_local{0};

  Array<Real> cache;

  Int release_{-1};
};

/* -------------------------------------------------------------------------- */
namespace internal {
  /* ------------------------------------------------------------------------ */
  template <class Array> class VectorPETScWrapped : public VectorPETSc {
  public:
    VectorPETScWrapped(Array && array) : array(array) {
      PETSc_call(VecCreateSeqWithArray, PETSC_COMM_SELF, 1, array.size(),
                 array.storage(), &x);
    }

    ~VectorPETScWrapped() { PETSc_call(VecDestroy, &x); }

  private:
    Array array;
  };

  /* ------------------------------------------------------------------------ */
  template <bool read_only> class VectorPETScLocal : public VectorPETSc {
  public:
    VectorPETScLocal(const Vec & g) : g(g) {
      PETSc_call(VecGetLocalVectorRead, g, x);
    }
    VectorPETScLocal(const VectorPETSc & g) : VectorPETScLocal(g.getVec()) {}
    ~VectorPETScLocal() {
      PETSc_call(VecRestoreLocalVectorRead, g, x);
      PETSc_call(VecDestroy, &x);
    }

    void printself(std::ostream & stream, int indent) const override {
      this->print(g, stream, indent);
    }

  private:
    const Vec & g;
  };

  template <> class VectorPETScLocal<false> : public VectorPETSc {
  public:
    VectorPETScLocal(Vec & g) : g(g) {
      PETSc_call(VecGetLocalVectorRead, g, x);
    }
    VectorPETScLocal(VectorPETSc & g) : VectorPETScLocal(g.getVec()) {}
    ~VectorPETScLocal() {
      PETSc_call(VecRestoreLocalVectorRead, g, x);
      PETSc_call(VecDestroy, &x);
    }

    void printself(std::ostream & stream, int indent) const override {
      this->print(g, stream, indent);
    }

  private:
    Vec & g;
  };

  /* ------------------------------------------------------------------------ */
  template <class Array>
  decltype(auto) make_petsc_wrapped_vector(Array && array) {
    return VectorPETScWrapped<Array>(std::forward<Array>(array));
  }

  template <
      typename V,
      std::enable_if_t<std::is_same<Vec, std::decay_t<V>>::value> * = nullptr>
  decltype(auto) make_petsc_local_vector(V && vec) {
    constexpr auto read_only = std::is_const<std::remove_reference_t<V>>::value;
    return VectorPETScLocal<read_only>(vec);
  }
} // namespace internal

} // namespace akantu

#include "vector_petsc_tmpl.hh"

#endif /* __AKANTU_SOLVER_VECTOR_PETSC_HH__ */
