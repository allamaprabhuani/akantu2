/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "aka_iterators.hh"
#include "aka_math.hh"
#include "aka_types.hh"
#include "communicator.hh"
/* -------------------------------------------------------------------------- */
#include <map>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_BBOX_HH_
#define AKANTU_AKA_BBOX_HH_

namespace akantu {

class BBox {
public:
  BBox() = default;

  BBox(Int spatial_dimension)
      : dim(spatial_dimension), lower_bounds(spatial_dimension),
        upper_bounds(spatial_dimension) {
    lower_bounds.fill(std::numeric_limits<Real>::max());
    upper_bounds.fill(-std::numeric_limits<Real>::max());
  }

  BBox(const BBox & other)
      : dim(other.dim), empty{false}, lower_bounds(other.lower_bounds),
        upper_bounds(other.upper_bounds) {}

  BBox & operator=(const BBox & other) {
    if (this != &other) {
      this->dim = other.dim;
      this->lower_bounds = other.lower_bounds;
      this->upper_bounds = other.upper_bounds;
      this->empty = other.empty;
    }
    return *this;
  }

  inline BBox & operator+=(const Vector<Real> & position) {
    AKANTU_DEBUG_ASSERT(
        this->dim == position.size(),
        "You are adding a point of a wrong dimension to the bounding box");

    this->empty = false;

    for (auto s : arange(dim)) {
      lower_bounds(s) = std::min(lower_bounds(s), position(s));
      upper_bounds(s) = std::max(upper_bounds(s), position(s));
    }
    return *this;
  }

  /* ------------------------------------------------------------------------ */
  inline bool intersects(const BBox & other,
                         const SpatialDirection & direction) const {
    AKANTU_DEBUG_ASSERT(
        this->dim == other.dim,
        "You are intersecting bounding boxes of different dimensions");
    return Math::intersects(lower_bounds(direction), upper_bounds(direction),
                            other.lower_bounds(direction),
                            other.upper_bounds(direction));
  }

  inline bool intersects(const BBox & other) const {
    if (this->empty or other.empty) {
      return false;
    }

    bool intersects_ = true;
    for (auto s : arange(this->dim)) {
      intersects_ &= this->intersects(other, SpatialDirection(s));
    }
    return intersects_;
  }

  /* ------------------------------------------------------------------------ */
  inline BBox intersection(const BBox & other) const {
    AKANTU_DEBUG_ASSERT(
        this->dim == other.dim,
        "You are intersecting bounding boxes of different dimensions");

    BBox intersection_(this->dim);
    intersection_.empty = not this->intersects(other);

    if (intersection_.empty) {
      return intersection_;
    }

    for (auto s : arange(this->dim)) {
      // is lower point in range ?
      bool point1 = Math::is_in_range(other.lower_bounds(s), lower_bounds(s),
                                      upper_bounds(s));

      // is upper point in range ?
      bool point2 = Math::is_in_range(other.upper_bounds(s), lower_bounds(s),
                                      upper_bounds(s));

      if (point1 and not point2) {
        // |-----------|         this (i)
        //       |-----------|   other(i)
        //       1           2
        intersection_.lower_bounds(s) = other.lower_bounds(s);
        intersection_.upper_bounds(s) = upper_bounds(s);
      } else if (point1 && point2) {
        // |-----------------|   this (i)
        //   |-----------|       other(i)
        //   1           2
        intersection_.lower_bounds(s) = other.lower_bounds(s);
        intersection_.upper_bounds(s) = other.upper_bounds(s);
      } else if (!point1 && point2) {
        //       |-----------|   this (i)
        // |-----------|         other(i)
        // 1           2
        intersection_.lower_bounds(s) = this->lower_bounds(s);
        intersection_.upper_bounds(s) = other.upper_bounds(s);
      } else {
        //   |-----------|       this (i)
        // |-----------------|   other(i)
        // 1                 2
        intersection_.lower_bounds(s) = this->lower_bounds(s);
        intersection_.upper_bounds(s) = this->upper_bounds(s);
      }
    }

    return intersection_;
  }

  /* ------------------------------------------------------------------------ */
  template <class Derived>
  inline bool contains(const Eigen::MatrixBase<Derived> & point) const {
    return (point.array() >= lower_bounds.array()).all() and
           (point.array() <= upper_bounds.array()).all();
  }

  /* ------------------------------------------------------------------------ */
  inline void reset() {
    lower_bounds.set(std::numeric_limits<Real>::max());
    upper_bounds.set(std::numeric_limits<Real>::lowest());
  }

  /* ------------------------------------------------------------------------ */
  template <class Derived>
  inline void getCenter(Eigen::MatrixBase<Derived> & center) {
    center = (upper_bounds + lower_bounds) / 2.;
  }

  /* ------------------------------------------------------------------------ */
  const Vector<Real> & getLowerBounds() const { return lower_bounds; }
  const Vector<Real> & getUpperBounds() const { return upper_bounds; }

  template <typename D>
  void setLowerBounds(const Eigen::MatrixBase<D> & lower_bounds) {
    this->lower_bounds = lower_bounds;
    this->empty = false;
  }
  template <typename D>
  void setUpperBounds(const Eigen::MatrixBase<D> & upper_bounds) {
    this->upper_bounds = upper_bounds;
    this->empty = false;
  }

  /* ------------------------------------------------------------------------ */
  inline Real size(const SpatialDirection & direction) const {
    return upper_bounds(direction) - lower_bounds(direction);
  }

  Vector<Real> size() const {
    Vector<Real> size_(dim);
    for (auto s : arange(this->dim)) {
      size_(s) = this->size(SpatialDirection(s));
    }
    return size_;
  }

  inline operator bool() const { return not empty; }

  /* ------------------------------------------------------------------------ */
  BBox allSum(const Communicator & communicator) const {
    Matrix<Real> reduce_bounds(dim, 2);

    reduce_bounds(0) = lower_bounds;
    reduce_bounds(1) = Real(-1.) * upper_bounds;

    communicator.allReduce(reduce_bounds, SynchronizerOperation::_min);

    BBox global(dim);
    global.lower_bounds = reduce_bounds(0);
    global.upper_bounds = Real(-1.) * reduce_bounds(1);
    global.empty = false;
    return global;
  }

  std::vector<BBox> allGather(const Communicator & communicator) const {
    auto prank = communicator.whoAmI();
    auto nb_proc = communicator.getNbProc();
    Array<Real> bboxes_data(nb_proc, dim * 2 + 1);

    auto * base = bboxes_data.data() + prank * (2 * dim + 1);
    MatrixProxy<Real> bounds(base, dim, 2);
    bounds(0) = lower_bounds;
    bounds(1) = upper_bounds;
    base[dim * 2] = empty ? 1. : 0.; // ugly trick

    communicator.allGather(bboxes_data);

    std::vector<BBox> bboxes;
    bboxes.reserve(nb_proc);

    for (auto p : arange(nb_proc)) {
      bboxes.emplace_back(dim);
      auto & bbox = bboxes.back();

      auto * base = bboxes_data.data() + p * (2 * dim + 1);
      MatrixProxy<Real> bounds(base, dim, 2);
      bbox.lower_bounds = bounds(0);
      bbox.upper_bounds = bounds(1);
      bbox.empty = (base[dim * 2] == 1.);
    }

    return bboxes;
  }

  std::map<Int, BBox> intersection(const BBox & other,
                                   const Communicator & communicator) const {
    // todo: change for a custom reduction algorithm
    auto other_bboxes = other.allGather(communicator);
    std::map<Int, BBox> intersections;
    for (const auto & bbox : enumerate(other_bboxes)) {
      auto && tmp = this->intersection(std::get<1>(bbox));
      if (tmp) {
        intersections[std::get<0>(bbox)] = tmp;
      }
    }
    return intersections;
  }

  void printself(std::ostream & stream) const {
    stream << "BBox[";
    if (not empty) {
      stream << lower_bounds << " - " << upper_bounds;
    }
    stream << "]";
  }

protected:
  Int dim{0};
  bool empty{true};
  Vector<Real> lower_bounds;
  Vector<Real> upper_bounds;
};

inline std::ostream & operator<<(std::ostream & stream, const BBox & bbox) {
  bbox.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* AKANTU_AKA_BBOX_HH_ */
