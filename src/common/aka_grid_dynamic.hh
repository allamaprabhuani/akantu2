/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_array.hh"
#include "aka_common.hh"
#include "aka_types.hh"
#include "mesh_accessor.hh"

#include <iostream>

/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_GRID_DYNAMIC_HH_
#define AKANTU_AKA_GRID_DYNAMIC_HH_

namespace akantu {

class Mesh;

template <typename T> class SpatialGrid {
public:
  explicit SpatialGrid(Int dimension)
      : dimension(dimension), spacing(dimension), center(dimension),
        lower(dimension), upper(dimension), empty_cell() {
    lower.fill(std::numeric_limits<Real>::max());
    upper.fill(-std::numeric_limits<Real>::max());
  }

  SpatialGrid(Int dimension, const Vector<Real> & spacing,
              const Vector<Real> & center)
      : SpatialGrid(dimension) {
    this->spacing = spacing;
    this->center = center;
  }

  virtual ~SpatialGrid() = default;

  class neighbor_cells_iterator;
  class cells_iterator;

  class CellID {
  public:
    CellID() = default;
    explicit CellID(Int dimention) : ids(dimention) {}
    void setID(Int dir, Int id) { ids(dir) = id; }
    Idx getID(Int dir) const { return ids(dir); }
    const Vector<Idx> & getIDs() const { return ids; }

    bool operator<(const CellID & id) const {
      return std::lexicographical_compare(ids.data(), ids.data() + ids.size(),
                                          id.ids.data(),
                                          id.ids.data() + id.ids.size());
    }

    bool operator==(const CellID & id) const {
      return std::equal(ids.data(), ids.data() + ids.size(), id.ids.data());
    }

    bool operator!=(const CellID & id) const { return !(operator==(id)); }

    class neighbor_cells_iterator {
    public:
      neighbor_cells_iterator(const CellID & cell_id, bool end)
          : cell_id(cell_id), position(cell_id.ids.size()) {
        position.fill(end ? 1 : -1);
        this->updateIt();
        if (end) {
          this->it++;
        }
      }

      neighbor_cells_iterator & operator++() {
        Int i = 0;
        for (; i < position.size() && position(i) == 1; ++i) {
        }

        if (i == position.size()) {
          ++it;
          return *this;
        }

        for (decltype(i) j = 0; j < i; ++j) {
          position(j) = -1;
        }
        position(i)++;
        updateIt();

        return *this;
      }

      neighbor_cells_iterator operator++(int) {
        neighbor_cells_iterator tmp(*this);
        operator++();
        return tmp;
      };

      bool operator==(const neighbor_cells_iterator & rhs) const {
        return cell_id == rhs.cell_id && it == rhs.it;
      };
      bool operator!=(const neighbor_cells_iterator & rhs) const {
        return !operator==(rhs);
      };

      CellID operator*() const {
        CellID cur_cell_id(cell_id);
        cur_cell_id.ids += position;
        return cur_cell_id;
      };

    private:
      void updateIt() {
        it = 0;
        for (Int i = 0; i < position.size(); ++i) {
          it = it * 3 + (position(i) + 1);
        }
      }

    private:
      /// central cell id
      const CellID & cell_id;
      // number representing the current neighbor in base 3;
      Int it{0};
      // current cell shift
      Vector<Idx> position;
    };

    class Neighbors {
    public:
      explicit Neighbors(const CellID & cell_id) : cell_id(cell_id) {}
      decltype(auto) begin() { return neighbor_cells_iterator(cell_id, false); }
      decltype(auto) end() { return neighbor_cells_iterator(cell_id, true); }

    private:
      const CellID & cell_id;
    };

    decltype(auto) neighbors() { return Neighbors(*this); }

  private:
    friend class cells_iterator;
    Vector<Idx> ids;
  };

  /* ------------------------------------------------------------------------ */
  class Cell {
  public:
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;

    Cell() : id(), data() {}

    explicit Cell(const CellID & cell_id) : id(cell_id), data() {}

    bool operator==(const Cell & cell) const { return id == cell.id; }
    bool operator!=(const Cell & cell) const { return id != cell.id; }

    Cell & add(const T & d) {
      data.push_back(d);
      return *this;
    }

    iterator begin() { return data.begin(); }
    const_iterator begin() const { return data.begin(); }

    iterator end() { return data.end(); }
    const_iterator end() const { return data.end(); }

  private:
    CellID id;
    std::vector<T> data;
  };

private:
  using cells_container = std::map<CellID, Cell>;

public:
  const Cell & getCell(const CellID & cell_id) const {
    auto it = cells.find(cell_id);
    if (it != cells.end()) {
      return it->second;
    }
    return empty_cell;
  }

  decltype(auto) beginCell(const CellID & cell_id) {
    auto it = cells.find(cell_id);
    if (it != cells.end()) {
      return it->second.begin();
    }
    return empty_cell.begin();
  }

  decltype(auto) endCell(const CellID & cell_id) {
    auto it = cells.find(cell_id);
    if (it != cells.end()) {
      return it->second.end();
    }
    return empty_cell.end();
  }

  decltype(auto) beginCell(const CellID & cell_id) const {
    auto it = cells.find(cell_id);
    if (it != cells.end()) {
      return it->second.begin();
    }
    return empty_cell.begin();
  }

  decltype(auto) endCell(const CellID & cell_id) const {
    auto it = cells.find(cell_id);
    if (it != cells.end()) {
      return it->second.end();
    }
    return empty_cell.end();
  }

  /* ------------------------------------------------------------------------ */
  class cells_iterator {
  public:
    explicit cells_iterator(typename std::map<CellID, Cell>::const_iterator it)
        : it(it) {}

    cells_iterator & operator++() {
      this->it++;
      return *this;
    }

    cells_iterator operator++(int /*unused*/) {
      cells_iterator tmp(*this);
      operator++();
      return tmp;
    };

    bool operator==(const cells_iterator & rhs) const { return it == rhs.it; };
    bool operator!=(const cells_iterator & rhs) const {
      return !operator==(rhs);
    };

    CellID operator*() const {
      CellID cur_cell_id(this->it->first);
      return cur_cell_id;
    };

  private:
    /// map iterator
    typename std::map<CellID, Cell>::const_iterator it;
  };

public:
  template <class vector_type>
  Cell & insert(const T & d, const vector_type & position) {
    auto && cell_id = getCellID(position);
    auto && it = cells.find(cell_id);
    if (it == cells.end()) {
      Cell cell(cell_id);
      auto & tmp = (cells[cell_id] = cell).add(d);

      auto posl =
          center.array() +
          cell_id.getIDs().array().template cast<Real>() * spacing.array();
      auto posu = posl + spacing.array();

      lower = lower.array().min(posl);
      upper = upper.array().max(posu);
      return tmp;
    }
    return it->second.add(d);
  }

  /* ------------------------------------------------------------------------ */
  inline decltype(auto) begin() const {
    auto begin = this->cells.begin();
    return cells_iterator(begin);
  }

  inline decltype(auto) end() const {
    auto end = this->cells.end();
    return cells_iterator(end);
  }

  template <class vector_type>
  CellID getCellID(const vector_type & position) const {
    CellID cell_id(dimension);
    for (Int i = 0; i < dimension; ++i) {
      cell_id.setID(i, getCellID(position(i), i));
    }
    return cell_id;
  }

  template <class vector_type>
  const Cell & getCell(const vector_type & position) const {
    return this->getCell(this->getCellID(position));
  }

  void printself(std::ostream & stream, int indent = 0) const {
    std::string space(indent, AKANTU_INDENT);

    std::streamsize prec = stream.precision();
    std::ios_base::fmtflags ff = stream.flags();

    stream.setf(std::ios_base::showbase);
    stream.precision(5);

    stream << space << "SpatialGrid<" << debug::demangle(typeid(T).name())
           << "> [" << std::endl;
    stream << space << " + dimension    : " << this->dimension << std::endl;
    stream << space << " + lower bounds : {";
    for (Int i = 0; i < lower.size(); ++i) {
      if (i != 0) {
        stream << ", ";
      }
      stream << lower(i);
    };
    stream << "}" << std::endl;
    stream << space << " + upper bounds : {";
    for (Int i = 0; i < upper.size(); ++i) {
      if (i != 0) {
        stream << ", ";
      }
      stream << upper(i);
    };
    stream << "}" << std::endl;
    stream << space << " + spacing : {";
    for (Int i = 0; i < spacing.size(); ++i) {
      if (i != 0) {
        stream << ", ";
      }
      stream << spacing(i);
    };
    stream << "}" << std::endl;
    stream << space << " + center : {";
    for (Int i = 0; i < center.size(); ++i) {
      if (i != 0) {
        stream << ", ";
      }
      stream << center(i);
    };
    stream << "}" << std::endl;

    stream << space << " + nb_cells     : " << this->cells.size() << "/";
    Vector<Real> dist(this->dimension);
    dist = upper;
    dist -= lower;
    for (Int i = 0; i < this->dimension; ++i) {
      dist(i) /= spacing(i);
    }
    auto nb_cells = std::ceil(dist(0));
    for (Int i = 1; i < this->dimension; ++i) {
      nb_cells *= std::ceil(dist(i));
    }
    stream << nb_cells << std::endl;
    stream << space << "]" << std::endl;

    stream.precision(prec);
    stream.flags(ff);
  }

  void saveAsMesh(Mesh & mesh) const;

private:
  /* --------------------------------------------------------------------------
   */
  inline decltype(auto) getCellID(Real position, Int direction) const {
    AKANTU_DEBUG_ASSERT(direction < center.size(), "The direction asked ("
                                                       << direction
                                                       << ") is out of range "
                                                       << center.size());
    Real dist_center = position - center(direction);
    Int id = std::floor(dist_center / spacing(direction));
    // if(dist_center < 0) id--;
    return id;
  }

  friend class GridSynchronizer;

public:
  AKANTU_GET_MACRO(LowerBounds, lower, const Vector<Real> &);
  AKANTU_GET_MACRO(UpperBounds, upper, const Vector<Real> &);
  AKANTU_GET_MACRO(Spacing, spacing, const Vector<Real> &);
  AKANTU_SET_MACRO(Spacing, spacing, Vector<Real> &);
  AKANTU_GET_MACRO(Center, center, const Vector<Real> &);
  AKANTU_SET_MACRO(Center, center, Vector<Real> &);

private:
  Int dimension{0};

  cells_container cells;

  Vector<Real> spacing;
  Vector<Real> center;

  Vector<Real> lower;
  Vector<Real> upper;

  Cell empty_cell;
};

/// standard output stream operator
template <typename T>
inline std::ostream & operator<<(std::ostream & stream,
                                 const SpatialGrid<T> & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "mesh.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T> void SpatialGrid<T>::saveAsMesh(Mesh & mesh) const {
  ElementType type = _not_defined;
  switch (dimension) {
  case 1:
    type = _segment_2;
    break;
  case 2:
    type = _quadrangle_4;
    break;
  case 3:
    type = _hexahedron_8;
    break;
  }

  MeshAccessor mesh_accessor(mesh);
  auto & connectivity = mesh_accessor.getConnectivity(type);
  auto & nodes = mesh_accessor.getNodes();
  auto & uint_data = mesh.getDataPointer<Int>("tag_1", type);

  Vector<Real> pos(dimension);

  Int global_id = 0;
  for (auto & cell_pair : cells) {
    auto cur_node = nodes.size();
    auto cur_elem = connectivity.size();
    const auto & cell_id = cell_pair.first;

    for (Int i = 0; i < dimension; ++i) {
      pos(i) = center(i) + cell_id.getID(i) * spacing(i);
    }
    nodes.push_back(pos);
    for (Int i = 0; i < dimension; ++i) {
      pos(i) += spacing(i);
    }
    nodes.push_back(pos);

    connectivity.push_back(cur_node);
    switch (dimension) {
    case 1:
      connectivity(cur_elem, 1) = cur_node + 1;
      break;
    case 2:
      pos(0) -= spacing(0);
      nodes.push_back(pos);
      pos(0) += spacing(0);
      pos(1) -= spacing(1);
      nodes.push_back(pos);
      connectivity(cur_elem, 1) = cur_node + 3;
      connectivity(cur_elem, 2) = cur_node + 1;
      connectivity(cur_elem, 3) = cur_node + 2;
      break;
    case 3:
      pos(1) -= spacing(1);
      pos(2) -= spacing(2);
      nodes.push_back(pos);
      pos(1) += spacing(1);
      nodes.push_back(pos);
      pos(0) -= spacing(0);
      nodes.push_back(pos);

      pos(1) -= spacing(1);
      pos(2) += spacing(2);
      nodes.push_back(pos);
      pos(0) += spacing(0);
      nodes.push_back(pos);
      pos(0) -= spacing(0);
      pos(1) += spacing(1);
      nodes.push_back(pos);

      connectivity(cur_elem, 1) = cur_node + 2;
      connectivity(cur_elem, 2) = cur_node + 3;
      connectivity(cur_elem, 3) = cur_node + 4;
      connectivity(cur_elem, 4) = cur_node + 5;
      connectivity(cur_elem, 5) = cur_node + 6;
      connectivity(cur_elem, 6) = cur_node + 1;
      connectivity(cur_elem, 7) = cur_node + 7;
      break;
    }
    uint_data.push_back(global_id);

    ++global_id;
  }
}

} // namespace akantu

#endif /* AKANTU_AKA_GRID_DYNAMIC_HH_ */
