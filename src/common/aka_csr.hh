/**
 * @file   aka_csr.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Apr 20 2011
 * @date last modification: Sun Oct 19 2014
 *
 * @brief  A compresed sparse row structure based on akantu Arrays
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "aka_array.hh"
#include "aka_common.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_CSR_HH__
#define __AKANTU_AKA_CSR_HH__

namespace akantu {

/**
 * This class  can be  used to  store the structure  of a  sparse matrix  or for
 * vectors with variable number of component per element
 *
 * @param nb_rows number of rows of a matrix or size of a vector.
 */
template <typename T> class CSR {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  explicit CSR(UInt nb_rows = 0)
      : nb_rows(nb_rows), rows_offsets(nb_rows + 1, 1, "rows_offsets"),
        rows(0, 1, "rows") {
    rows_offsets.clear();
  };

  virtual ~CSR() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// does nothing
  inline void beginInsertions(){};

  /// insert a new entry val in row row
  inline UInt insertInRow(UInt row, const T & val) {
    UInt pos = rows_offsets(row)++;
    rows(pos) = val;
    return pos;
  }

  /// access an element of the matrix
  inline const T & operator()(UInt row, UInt col) const {
    AKANTU_DEBUG_ASSERT(rows_offsets(row + 1) - rows_offsets(row) > col,
                        "This element is not present in this CSR");
    return rows(rows_offsets(row) + col);
  }

  /// access an element of the matrix
  inline T & operator()(UInt row, UInt col) {
    AKANTU_DEBUG_ASSERT(rows_offsets(row + 1) - rows_offsets(row) > col,
                        "This element is not present in this CSR");
    return rows(rows_offsets(row) + col);
  }

  inline void endInsertions() {
    for (UInt i = nb_rows; i > 0; --i)
      rows_offsets(i) = rows_offsets(i - 1);
    rows_offsets(0) = 0;
  }

  inline void countToCSR() {
    for (UInt i = 1; i < nb_rows; ++i)
      rows_offsets(i) += rows_offsets(i - 1);
    for (UInt i = nb_rows; i >= 1; --i)
      rows_offsets(i) = rows_offsets(i - 1);
    rows_offsets(0) = 0;
  }

  inline void clearRows() {
    rows_offsets.clear();
    rows.resize(0);
  };

  inline void resizeRows(UInt nb_rows) {
    this->nb_rows = nb_rows;
    rows_offsets.resize(nb_rows + 1);
    rows_offsets.clear();
  }

  inline void resizeCols() { rows.resize(rows_offsets(nb_rows)); }

  inline void copy(Array<UInt> & offsets, Array<T> & values) {
    offsets.copy(rows_offsets);
    values.copy(rows);
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// returns the number of rows
  inline UInt getNbRows() const { return rows_offsets.size() - 1; };

  /// returns the number of non-empty columns in a given row
  inline UInt getNbCols(UInt row) const {
    return rows_offsets(row + 1) - rows_offsets(row);
  };

  /// returns the offset (start of columns) for a given row
  inline UInt & rowOffset(UInt row) { return rows_offsets(row); };

  /// iterator on a row
  template <class R>
  class iterator_internal
      : public std::iterator<std::bidirectional_iterator_tag, R> {
  public:
    using _parent = std::iterator<std::bidirectional_iterator_tag, R>;
    using pointer = typename _parent::pointer;
    using reference = typename _parent::reference;

    explicit iterator_internal(pointer x = nullptr) : pos(x){};
    iterator_internal(const iterator_internal & it) : pos(it.pos){};

    iterator_internal & operator++() {
      ++pos;
      return *this;
    };
    iterator_internal operator++(int) {
      iterator tmp(*this);
      operator++();
      return tmp;
    };

    iterator_internal & operator--() {
      --pos;
      return *this;
    };
    iterator_internal operator--(int) {
      iterator_internal tmp(*this);
      operator--();
      return tmp;
    };

    bool operator==(const iterator_internal & rhs) { return pos == rhs.pos; };
    bool operator!=(const iterator_internal & rhs) { return pos != rhs.pos; };
    reference operator*() { return *pos; };
    pointer operator->() const { return pos; };

  private:
    pointer pos;
  };

  using iterator = iterator_internal<T>;
  using const_iterator = iterator_internal<const T>;

  inline iterator begin(UInt row) {
    return iterator(rows.storage() + rows_offsets(row));
  };
  inline iterator end(UInt row) {
    return iterator(rows.storage() + rows_offsets(row + 1));
  };

  inline const_iterator begin(UInt row) const {
    return const_iterator(rows.storage() + rows_offsets(row));
  };
  inline const_iterator end(UInt row) const {
    return const_iterator(rows.storage() + rows_offsets(row + 1));
  };

  inline iterator rbegin(UInt row) {
    return iterator(rows.storage() + rows_offsets(row + 1) - 1);
  };
  inline iterator rend(UInt row) {
    return iterator(rows.storage() + rows_offsets(row) - 1);
  };

  inline const Array<UInt> & getRowsOffset() const { return rows_offsets; };
  inline const Array<T> & getRows() const { return rows; };
  inline Array<T> & getRows() { return rows; };
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  UInt nb_rows;

  /// array of size nb_rows containing the offset where the values are stored in
  Array<UInt> rows_offsets;

  /// compressed row values, values of row[i] are stored between rows_offsets[i]
  /// and rows_offsets[i+1]
  Array<T> rows;
};

/* -------------------------------------------------------------------------- */
/* Data CSR                                                                   */
/* -------------------------------------------------------------------------- */

/**
 * Inherits from  CSR<UInt> and  can contain information  such as  matrix values
 * where the mother class would be a CSR structure for row and cols
 *
 * @return nb_rows
 */
template <class T> class DataCSR : public CSR<UInt> {
public:
  DataCSR(UInt nb_rows = 0) : CSR<UInt>(nb_rows), data(0, 1){};

  inline void resizeCols() {
    CSR<UInt>::resizeCols();
    data.resize(rows_offsets(nb_rows));
  }

  inline const Array<T> & getData() const { return data; };

private:
  Array<T> data;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "aka_csr_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const CSR & _this)
// {
//   _this.printself(stream);
//   return stream;
// }

} // akantu

#endif /* __AKANTU_AKA_CSR_HH__ */
