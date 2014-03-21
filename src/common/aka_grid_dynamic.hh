/**
 * @file   aka_grid_dynamic.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Feb 15 16:55:32 2013
 *
 * @brief  Grid that is auto balanced
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
#include "aka_common.hh"
#include "aka_array.hh"
#include "aka_types.hh"

#include <iostream>

/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_GRID_DYNAMIC_HH__
#define __AKANTU_AKA_GRID_DYNAMIC_HH__

__BEGIN_AKANTU__

class Mesh;

template<typename T>
class SpatialGrid {
public:
  SpatialGrid(UInt dimension) : dimension(dimension),
                                spacing(dimension),
                                center(dimension),
                                lower(dimension),
                                upper(dimension),
                                empty_cell() {}

  SpatialGrid(UInt dimension,
              const Vector<Real> & spacing,
              const Vector<Real> & center) : dimension(dimension),
                                                    spacing(spacing),
                                                    center(center),
                                                    lower(dimension),
                                                    upper(dimension),
                                                    empty_cell() {
    for (UInt i = 0; i < dimension; ++i) {
      lower(i) =   std::numeric_limits<Real>::max();
      upper(i) = - std::numeric_limits<Real>::max();
    }
  }

  virtual ~SpatialGrid() {};

  class neighbor_cells_iterator;

  class CellID {
  public:
    CellID() : ids() {}
    CellID(UInt dimention) : ids(dimention) {}
    void setID(UInt dir, Int id) { ids(dir) = id; }
    Int getID(UInt dir) const { return ids(dir); }

    bool operator<(const CellID & id) const {
      return std::lexicographical_compare(ids.storage(), ids.storage() + ids.size(),
                                          id.ids.storage(), id.ids.storage() + id.ids.size());
    }

    bool operator==(const CellID & id) const {
      return std::equal(ids.storage(), ids.storage() + ids.size(),
                        id.ids.storage());
    }

    bool operator!=(const CellID & id) const {
      return !(operator==(id));
    }

  private:
    friend class neighbor_cells_iterator;
    Vector<Int> ids;
  };

  /* -------------------------------------------------------------------------- */
  class Cell {
  public:
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;

    Cell() : id(), data() { }

    Cell(const CellID &cell_id) : id(cell_id), data() { }

    bool operator==(const Cell & cell) const { return id == cell.id; }
    bool operator!=(const Cell & cell) const { return id != cell.id; }

    Cell & add(const T & d) { data.push_back(d); return *this; }

    iterator begin() { return data.begin(); }
    const_iterator begin() const { return data.begin(); }

    iterator end() { return data.end(); }
    const_iterator end() const { return data.end(); }

#if not defined(AKANTU_NDEBUG)
    Cell & add(const T & d, const Vector<Real> & pos) {
      data.push_back(d); positions.push_back(pos); return *this;
    }
    typedef typename std::vector< Vector<Real> >::const_iterator position_iterator;
    position_iterator begin_pos() const { return positions.begin(); }
    position_iterator end_pos() const { return positions.end(); }
#endif

  private:
    CellID id;
    std::vector<T> data;
#if not defined(AKANTU_NDEBUG)
    std::vector< Vector<Real> > positions;
#endif
  };

private:
  typedef std::map<CellID, Cell> cells_container;

public:

  const Cell & getCell(const CellID & cell_id) const {
    typename cells_container::const_iterator it = cells.find(cell_id);
    if(it != cells.end()) return it->second;
    else return empty_cell;
  }

  typename Cell::iterator beginCell(const CellID & cell_id) {
    typename cells_container::iterator it = cells.find(cell_id);
    if(it != cells.end()) return it->second.begin();
    else return empty_cell.begin();
  }

  typename Cell::iterator endCell(const CellID & cell_id) {
    typename cells_container::iterator it = cells.find(cell_id);
    if(it != cells.end()) return it->second.end();
    else return empty_cell.end();
  }

  typename Cell::const_iterator beginCell(const CellID & cell_id) const {
    typename cells_container::const_iterator it = cells.find(cell_id);
    if(it != cells.end()) return it->second.begin();
    else return empty_cell.begin();
  }

  typename Cell::const_iterator endCell(const CellID & cell_id) const {
    typename cells_container::const_iterator it = cells.find(cell_id);
    if(it != cells.end()) return it->second.end();
    else return empty_cell.end();
  }


  class neighbor_cells_iterator : private std::iterator<std::forward_iterator_tag, UInt> {
  public:
    neighbor_cells_iterator(const CellID & cell_id, bool end) :
      cell_id(cell_id),
      position(cell_id.ids.size(), end ? 1 : -1) {

      this->updateIt();
      if(end) this->it++;
    }

    neighbor_cells_iterator& operator++() {
      UInt i = 0;
      for (; i < position.size() && position(i) == 1; ++i);

      if(i == position.size()) ++it;
      else {
        for (UInt j = 0; j < i; ++j) position(j) = -1;
        position(i)++;
        updateIt();
      }

      return *this;
    }

    neighbor_cells_iterator operator++(int) { neighbor_cells_iterator tmp(*this); operator++(); return tmp; };

    bool operator==(const neighbor_cells_iterator& rhs) const { return cell_id == rhs.cell_id && it == rhs.it; };
    bool operator!=(const neighbor_cells_iterator& rhs) const { return ! operator==(rhs); };

    CellID operator*() const {
      CellID cur_cell_id(cell_id);
      cur_cell_id.ids += position;
      return cur_cell_id;
    };

  private:
    void updateIt() {
      it = 0;
      for (UInt i = 0; i < position.size(); ++i) it = it * 3 + (position(i) + 1);
    }

  private:
    /// central cell id
    const CellID & cell_id;

    // number representing the current neighbor in base 3;
    UInt it;

    Vector<Int> position;
  };

public:
  template<class vector_type>
  Cell & insert(const T & d, const vector_type & position) {
    CellID cell_id = getCellID(position);

    typename cells_container::iterator it = cells.find(cell_id);

    if(it == cells.end()) {
      Cell cell(cell_id);
#if defined(AKANTU_NDEBUG)
      Cell & tmp = (cells[cell_id] = cell).add(d);
#else
      Cell & tmp = (cells[cell_id] = cell).add(d, position);
#endif

      for (UInt i = 0; i < dimension; ++i) {
        Real posl = center(i) + cell_id.getID(i) * spacing(i);
        Real posu = posl + spacing(i);
        if(posl < lower(i)) lower(i) = posl;
        if(posu > upper(i)) upper(i) = posu;
      }
      return tmp;
    } else {
#if defined(AKANTU_NDEBUG)
      return it->second.add(d);
#else
      return it->second.add(d, position);
#endif
    }
  }


  inline neighbor_cells_iterator beginNeighborCells(const CellID & cell_id) const {
    return neighbor_cells_iterator(cell_id, false);
  }

  inline neighbor_cells_iterator endNeighborCells(const CellID & cell_id) const {
    return neighbor_cells_iterator(cell_id, true);
  }


  template<class vector_type>
  CellID getCellID(const vector_type & position) const {
    CellID cell_id(dimension);
    for (UInt i = 0; i < dimension; ++i) {
      cell_id.setID(i, getCellID(position(i), i));
    }
    return cell_id;
  }


  void printself(std::ostream & stream, int indent = 0) const {
    std::string space;
    for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

    std::streamsize prec        = stream.precision();
    std::ios_base::fmtflags ff  = stream.flags();

    stream.setf (std::ios_base::showbase);
    stream.precision(5);

    stream << space << "SpatialGrid<" << debug::demangle(typeid(T).name()) << "> [" << std::endl;
    stream << space << " + dimension    : " << this->dimension << std::endl;
    stream << space << " + lower bounds : {";
    for (UInt i = 0; i < lower.size(); ++i) { if(i != 0) stream << ", "; stream << lower(i); };
    stream << "}" << std::endl;
    stream << space << " + upper bounds : {";
    for (UInt i = 0; i < upper.size(); ++i) { if(i != 0) stream << ", "; stream << upper(i); };
    stream << "}" << std::endl;
    stream << space << " + spacing : {";
    for (UInt i = 0; i < spacing.size(); ++i) { if(i != 0) stream << ", "; stream << spacing(i); };
    stream << "}" << std::endl;
    stream << space << " + center : {";
    for (UInt i = 0; i < center.size(); ++i) { if(i != 0) stream << ", "; stream << center(i); };
    stream << "}" << std::endl;

    stream << space << " + nb_cells     : " << this->cells.size() << "/";
    Vector<Real> dist(this->dimension);
    dist = upper;
    dist -= lower;
    for (UInt i = 0; i < this->dimension; ++i) {
      dist(i) /= spacing(i); 
    }
    UInt nb_cells = std::ceil(dist(0));
    for (UInt i = 1; i < this->dimension; ++i) {
      nb_cells *= std::ceil(dist(i));
    }
    stream << nb_cells << std::endl;
    stream << space << "]" << std::endl;

    stream.precision(prec);
    stream.flags(ff);
  }

  void saveAsMesh(Mesh & mesh) const;

private:
  /* -------------------------------------------------------------------------- */
  inline UInt getCellID(Real position, UInt direction) const {
    AKANTU_DEBUG_ASSERT(direction < center.size(), "The direction asked ("
                        << direction << ") is out of range " << center.size());
    Real dist_center = position - center(direction);
    Int id = std::floor(dist_center / spacing(direction));
    //if(dist_center < 0) id--;
    return id;
  }

  friend class GridSynchronizer;
public:
  AKANTU_GET_MACRO(LowerBounds, lower, const Vector<Real> &);
  AKANTU_GET_MACRO(UpperBounds, upper, const Vector<Real> &);
  AKANTU_GET_MACRO(Spacing, spacing, const Vector<Real> &);

protected:
  UInt dimension;

  cells_container cells;

  Vector<Real> spacing;
  Vector<Real> center;

  Vector<Real> lower;
  Vector<Real> upper;

  Cell empty_cell;
};

/// standard output stream operator
template<typename T>
inline std::ostream & operator <<(std::ostream & stream, const SpatialGrid<T> & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__
#include "mesh.hh"
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<typename T>
void SpatialGrid<T>::saveAsMesh(Mesh & mesh) const {
  Array<Real> & nodes = const_cast<Array<Real> &>(mesh.getNodes());

  ElementType type;
  switch(dimension) {
  case 1: type = _segment_2; break;
  case 2: type = _quadrangle_4; break;
  case 3: type = _hexahedron_8; break;
  }

  mesh.addConnectivityType(type);
  Array<UInt> & connectivity = const_cast<Array<UInt> &>(mesh.getConnectivity(type));
  Array<UInt> & uint_data = *mesh.getDataPointer<UInt>("tag_1", type);

  typename cells_container::const_iterator it  = cells.begin();
  typename cells_container::const_iterator end = cells.end();

  Vector<Real> pos(dimension);

  UInt global_id = 0;
  for (;it != end; ++it, ++global_id) {
    UInt cur_node = nodes.getSize();
    UInt cur_elem = connectivity.getSize();
    const CellID & cell_id = it->first;

    for (UInt i = 0; i < dimension; ++i)  pos(i) = center(i) + cell_id.getID(i) * spacing(i);
    nodes.push_back(pos);
    for (UInt i = 0; i < dimension; ++i)  pos(i) += spacing(i);
    nodes.push_back(pos);

    connectivity.push_back(cur_node);
    switch(dimension) {
    case 1:
      connectivity(cur_elem, 1) = cur_node + 1;
      break;
    case 2:
      pos(0) -= spacing(0); nodes.push_back(pos);
      pos(0) += spacing(0); pos(1) -= spacing(1); nodes.push_back(pos);
      connectivity(cur_elem, 1) = cur_node + 3;
      connectivity(cur_elem, 2) = cur_node + 1;
      connectivity(cur_elem, 3) = cur_node + 2;
      break;
    case 3:
                            pos(1) -= spacing(1); pos(2) -= spacing(2); nodes.push_back(pos);
                            pos(1) += spacing(1);                       nodes.push_back(pos);
      pos(0) -= spacing(0);                                             nodes.push_back(pos);

                            pos(1) -= spacing(1); pos(2) += spacing(2); nodes.push_back(pos);
      pos(0) += spacing(0);                                             nodes.push_back(pos);
      pos(0) -= spacing(0); pos(1) += spacing(1);                       nodes.push_back(pos);


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
  }

#if not defined(AKANTU_NDEBUG)
  mesh.addConnectivityType(_point_1);
  Array<UInt> & connectivity_pos = const_cast<Array<UInt> &>(mesh.getConnectivity(_point_1));
  Array<UInt> & uint_data_pos = *mesh.getDataPointer<UInt>( "tag_1", _point_1);
  Array<UInt> & uint_data_pos_ghost = *mesh.getDataPointer<UInt>("tag_0", _point_1);

  it  = cells.begin();
  global_id = 0;
  for (;it != end; ++it, ++global_id) {
    typename Cell::position_iterator cell_it  = it->second.begin_pos();
    typename Cell::const_iterator    cell_it_cont  = it->second.begin();
    typename Cell::position_iterator cell_end = it->second.end_pos();
    for (;cell_it != cell_end; ++cell_it, ++cell_it_cont) {
      nodes.push_back(*cell_it);
      connectivity_pos.push_back(nodes.getSize()-1);
      uint_data_pos.push_back(global_id);
      uint_data_pos_ghost.push_back(cell_it_cont->ghost_type==_ghost);
    }
  }
#endif
}

__END_AKANTU__



#endif /* __AKANTU_AKA_GRID_DYNAMIC_HH__ */
