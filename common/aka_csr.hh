/**
 * @file   aka_csr.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Apr 18 13:39:52 2011
 *
 * @brief  A compresed sparse row structure based on akantu Vectors
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

#ifndef __AKANTU_AKA_CSR_HH__
#define __AKANTU_AKA_CSR_HH__

__BEGIN_AKANTU__

template <typename T>
class CSR {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  CSR(UInt nb_rows = 0) : nb_rows(nb_rows), rows_offsets(nb_rows + 1, 1), rows(0,1) { rows_offsets.clear(); };
  virtual ~CSR() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  inline void beginInsertions() {};

  inline UInt insertInRow(UInt row, T & val) { 
    UInt pos = rows_offsets(row)++;
    rows(pos) = val;
    return pos;
  }

  inline const T & operator()(UInt row, UInt col) const {
    return rows(rows_offsets(row) + col);
  }

  inline T & operator()(UInt row, UInt col) {
    return rows(rows_offsets(row) + col);
  }

  inline void endInsertions() {
    for (UInt i = nb_rows; i > 0; --i) rows_offsets(i) = rows_offsets(i-1);
    rows_offsets(0) = 0;
  }

  inline void countToCSR() {
    for (UInt i = 1; i < nb_rows; ++i) rows_offsets(i) += rows_offsets(i-1);
    for (UInt i = nb_rows; i >= 1; --i) rows_offsets(i) = rows_offsets(i-1);
    rows_offsets(0) = 0;
  }

  inline void clearRows() { rows_offsets.clear(); rows.resize(0); };

  inline void resizeRows(UInt nb_rows) { 
    this->nb_rows = nb_rows;
    rows_offsets.resize(nb_rows + 1);
  }

  inline void resizeCols() { 
    rows.resize(rows_offsets(nb_rows));
  }

  inline void copy(Vector<UInt> & offsets, Vector<T> & values) {
    offsets.copy(rows_offsets);
    values.copy(rows);
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  inline UInt getNbRows() { return rows_offsets.getSize() - 1; };

  inline UInt getNbCols(UInt row) { return rows_offsets(row + 1) - rows_offsets(row); };

  inline UInt & rowOffset(UInt row) { return rows_offsets(row); };

  /// iterator on a row
  class iterator : public std::iterator<std::bidirectional_iterator_tag, T> {
  public:
    typedef std::iterator<std::input_iterator_tag, T> _parent;
    typedef typename _parent::pointer pointer;
    typedef typename _parent::reference reference;
    iterator(pointer x = NULL) : pos(x) {};
    iterator(const iterator & it) : pos(it.pos) {};

    iterator& operator++() { ++pos; return *this; };
    iterator operator++(int) { iterator tmp(*this); operator++(); return tmp; };

    iterator& operator--() { --pos; return *this; };
    iterator operator--(int) { iterator tmp(*this); operator--(); return tmp; };

    bool operator==(const iterator& rhs) { return pos == rhs.pos; };
    bool operator!=(const iterator& rhs) { return pos != rhs.pos; };
    reference operator*() { return *pos; };
    pointer operator->() const { return pos; };
  private:
    pointer pos;
  };

  inline iterator begin(UInt row) { return iterator(rows.values + rows_offsets(row)); };
  inline iterator end(UInt row) { return iterator(rows.values + rows_offsets(row+1)); };

  inline iterator rbegin(UInt row) { return iterator(rows.values + rows_offsets(row+1) - 1); };
  inline iterator rend(UInt row) { return iterator(rows.values + rows_offsets(row) - 1); };

  inline const Vector<UInt> & getRowsOffset() const { return rows_offsets; };
  inline const Vector<T> & getRows() const { return rows; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  UInt nb_rows;

  /// array of size nb_rows containing the offset where the values are stored in
  Vector<UInt> rows_offsets;

  /// compressed row values, values of row[i] are stored between rows_offsets[i]
  /// and rows_offsets[i+1]
  Vector<T> rows;
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


__END_AKANTU__

#endif /* __AKANTU_AKA_CSR_HH__ */
