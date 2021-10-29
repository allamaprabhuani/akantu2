/**
 * @file   aka_array_tmpl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jul 15 2010
 * @date last modification: Fri Feb 26 2021
 *
 * @brief  Inline functions of the classes Array<T> and ArrayBase
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
/* Inline Functions Array<T>                                                  */
/* -------------------------------------------------------------------------- */
#include "aka_array.hh" // NOLINT
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */
//#ifndef __AKANTU_AKA_ARRAY_TMPL_HH__
//#define __AKANTU_AKA_ARRAY_TMPL_HH__
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace debug {
  struct ArrayException : public Exception {};
} // namespace debug

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
ArrayDataLayer<T, allocation_trait>::ArrayDataLayer(Int size, Int nb_component,
                                                    const ID & id)
    : ArrayBase(id) {
  allocate(size, nb_component);
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
ArrayDataLayer<T, allocation_trait>::ArrayDataLayer(Int size, Int nb_component,
                                                    const_reference value,
                                                    const ID & id)
    : ArrayBase(id) {
  allocate(size, nb_component, value);
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
ArrayDataLayer<T, allocation_trait>::ArrayDataLayer(const ArrayDataLayer & vect,
                                                    const ID & id)
    : ArrayBase(vect, id) {
  this->data_storage = vect.data_storage;
  this->size_ = vect.size_;
  this->nb_component = vect.nb_component;
  this->values = this->data_storage.data();
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
ArrayDataLayer<T, allocation_trait>::ArrayDataLayer(
    const std::vector<value_type> & vect) {
  this->data_storage = vect;
  this->size_ = vect.size();
  this->nb_component = 1;
  this->values = this->data_storage.data();
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
ArrayDataLayer<T, allocation_trait> &
ArrayDataLayer<T, allocation_trait>::operator=(const ArrayDataLayer & other) {
  if (this != &other) {
    this->data_storage = other.data_storage;
    this->nb_component = other.nb_component;
    this->size_ = other.size_;
    this->values = this->data_storage.data();
  }
  return *this;
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
ArrayDataLayer<T, allocation_trait>::ArrayDataLayer(ArrayDataLayer && other) =
    default;

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
ArrayDataLayer<T, allocation_trait> &
ArrayDataLayer<T, allocation_trait>::operator=(ArrayDataLayer && other) =
    default;

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
void ArrayDataLayer<T, allocation_trait>::allocate(Int new_size,
                                                   Int nb_component) {
  this->nb_component = nb_component;
  this->resize(new_size);
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
void ArrayDataLayer<T, allocation_trait>::allocate(Int new_size,
                                                   Int nb_component,
                                                   const T & val) {
  this->nb_component = nb_component;
  this->resize(new_size, val);
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
void ArrayDataLayer<T, allocation_trait>::resize(Int new_size) {
  this->data_storage.resize(new_size * this->nb_component);
  this->values = this->data_storage.data();
  this->size_ = new_size;
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
void ArrayDataLayer<T, allocation_trait>::resize(Int new_size,
                                                 const T & value) {
  this->data_storage.resize(new_size * this->nb_component, value);
  this->values = this->data_storage.data();
  this->size_ = new_size;
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
void ArrayDataLayer<T, allocation_trait>::reserve(Int size, Int new_size) {
  if (new_size != -1) {
    this->data_storage.resize(new_size * this->nb_component);
  }

  this->data_storage.reserve(size * this->nb_component);
  this->values = this->data_storage.data();
}

/* -------------------------------------------------------------------------- */
/**
 * append a tuple to the array with the value value for all components
 * @param value the new last tuple or the array will contain nb_component copies
 * of value
 */
template <typename T, ArrayAllocationType allocation_trait>
inline void ArrayDataLayer<T, allocation_trait>::push_back(const T & value) {
  this->data_storage.push_back(value);
  this->values = this->data_storage.data();
  this->size_ += 1;
}

/* -------------------------------------------------------------------------- */
/**
 * append a matrix or a vector to the array
 * @param new_elem a reference to a Matrix<T> or Vector<T> */
template <typename T, ArrayAllocationType allocation_trait>
template <typename Derived>
inline void ArrayDataLayer<T, allocation_trait>::push_back(
    const Eigen::MatrixBase<Derived> & new_elem) {
  AKANTU_DEBUG_ASSERT(
      nb_component == new_elem.size(),
      "The vector("
          << new_elem.size()
          << ") as not a size compatible with the Array (nb_component="
          << nb_component << ").");
  for (Idx i = 0; i < new_elem.size(); ++i) {
    this->data_storage.push_back(new_elem.data()[i]);
  }
  this->values = this->data_storage.data();
  this->size_ += 1;
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
inline Int ArrayDataLayer<T, allocation_trait>::getAllocatedSize() const {
  return this->data_storage.capacity() / this->nb_component;
}

/* -------------------------------------------------------------------------- */
template <typename T, ArrayAllocationType allocation_trait>
inline Int ArrayDataLayer<T, allocation_trait>::getMemorySize() const {
  return this->data_storage.capacity() * sizeof(T);
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <typename T>
class ArrayDataLayer<T, ArrayAllocationType::_pod> : public ArrayBase {
public:
  using value_type = T;
  using reference = value_type &;
  using pointer_type = value_type *;
  using const_reference = const value_type &;

public:
  ~ArrayDataLayer() override { deallocate(); }

  /// Allocation of a new vector
  ArrayDataLayer(Int size = 0, Int nb_component = 1, const ID & id = "")
      : ArrayBase(id) {
    allocate(size, nb_component);
  }

  /// Allocation of a new vector with a default value
  ArrayDataLayer(Int size, Int nb_component, const_reference value,
                 const ID & id = "")
      : ArrayBase(id) {
    allocate(size, nb_component, value);
  }

  /// Copy constructor (deep copy)
  ArrayDataLayer(const ArrayDataLayer & vect, const ID & id = "")
      : ArrayBase(vect, id) {
    allocate(vect.size(), vect.getNbComponent());
    std::copy_n(vect.data(), this->size_ * this->nb_component, values);
  }

  /// Copy constructor (deep copy)
  explicit ArrayDataLayer(const std::vector<value_type> & vect) {
    allocate(vect.size(), 1);
    std::copy_n(vect.data(), this->size_ * this->nb_component, values);
  }

  // copy operator
  inline ArrayDataLayer & operator=(const ArrayDataLayer & other) {
    if (this != &other) {
      allocate(other.size(), other.getNbComponent());
      std::copy_n(other.data(), this->size_ * this->nb_component, values);
    }
    return *this;
  }

  // move constructor
  inline ArrayDataLayer(ArrayDataLayer && other) noexcept = default;

  // move assign
  inline ArrayDataLayer & operator=(ArrayDataLayer && other) noexcept = default;

protected:
  // deallocate the memory
  virtual void deallocate() {
    // NOLINTNEXTLINE(cppcoreguidelines-owning-memory,
    // cppcoreguidelines-no-malloc)
    free(this->values);
  }

  // allocate the memory
  virtual inline void allocate(Int size, Int nb_component) {
    if (size != 0) { // malloc can return a non NULL pointer in case size is 0
      this->values = static_cast<T *>(                   // NOLINT
          std::malloc(nb_component * size * sizeof(T))); // NOLINT
    }

    if (this->values == nullptr and size != 0) {
      throw std::bad_alloc();
    }
    this->nb_component = nb_component;
    this->allocated_size = this->size_ = size;
  }

  // allocate and initialize the memory
  virtual inline void allocate(Int size, Int nb_component, const T & value) {
    allocate(size, nb_component);
    std::fill_n(values, size * nb_component, value);
  }

public:
  /// append a tuple of size nb_component containing value
  inline void push_back(const_reference value) {
    resize(this->size_ + 1, value);
  }

  /// append a Vector or a Matrix
  template <typename Derived>
  inline void push_back(const Eigen::MatrixBase<Derived> & new_elem) {
    AKANTU_DEBUG_ASSERT(
        nb_component == new_elem.size(),
        "The vector("
            << new_elem.size()
            << ") as not a size compatible with the Array (nb_component="
            << nb_component << ").");
    this->resize(this->size_ + 1);
    std::copy_n(new_elem.derived().data(), new_elem.size(),
                values + this->nb_component * (this->size_ - 1));
  }

  /// changes the allocated size but not the size
  virtual void reserve(Int size, Int new_size = Int(-1)) {
    auto tmp_size = this->size_;
    if (new_size != Int(-1))
      tmp_size = new_size;
  }
  this->resize(size);
  this->size_ = std::min(this->size_, tmp_size);
}

/// change the size of the Array
virtual void
resize(Int size) {
  if (size * this->nb_component == 0) {
    free(values); // NOLINT: cppcoreguidelines-no-malloc
    values = nullptr;
    this->allocated_size = 0;
  } else {
    if (this->values == nullptr) {
      this->allocate(size, this->nb_component);
      return;
    }

    Int diff = size - allocated_size;
    Int size_to_allocate = (std::abs(diff) > AKANTU_MIN_ALLOCATION) ? size
                           : (diff > 0) ? allocated_size + AKANTU_MIN_ALLOCATION
                                        : allocated_size;

    if (size_to_allocate ==
        allocated_size) { // otherwhy the reserve + push_back might fail...
      this->size_ = size;
      return;
    }

    auto * tmp_ptr = reinterpret_cast<T *>( // NOLINT
        realloc(this->values,
                size_to_allocate * this->nb_component * sizeof(T)));

    if (tmp_ptr == nullptr) {
      throw std::bad_alloc();
    }

    this->values = tmp_ptr;
    this->allocated_size = size_to_allocate;
  }

  this->size_ = size;
}

/// change the size of the Array and initialize the values
virtual void resize(Int size, const T & val) {
  Int tmp_size = this->size_;
  this->resize(size);
  if (size > tmp_size) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    std::fill_n(values + this->nb_component * tmp_size,
                (size - tmp_size) * this->nb_component, val);
  }
}

/// get the amount of space allocated in bytes
inline UInt getMemorySize() const final {
  return this->allocated_size * this->nb_component * sizeof(T);
}

/// Get the real size allocated in memory
inline Int getAllocatedSize() const { return this->allocated_size; }

/// give the address of the memory allocated for this vector
[[deprecated("use data instead to be stl compatible")]] T * storage() const {
  return values;
};
T * data() const { return values; };

protected:
/// allocation type agnostic  data access
T * values{nullptr};

Int allocated_size{0};
}; // namespace akantu

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline auto Array<T, is_scal>::operator()(Int i, Int j) -> reference {
  AKANTU_DEBUG_ASSERT(this->size_ > 0,
                      "The array \"" << this->id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < this->size_) && (j < this->nb_component),
                      "The value at position ["
                          << i << "," << j << "] is out of range in array \""
                          << this->id << "\"");
  return this->values[i * this->nb_component + j];
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline auto Array<T, is_scal>::operator()(Int i, Int j) const
    -> const_reference {
  AKANTU_DEBUG_ASSERT(this->size_ > 0,
                      "The array \"" << this->id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < this->size_) && (j < this->nb_component),
                      "The value at position ["
                          << i << "," << j << "] is out of range in array \""
                          << this->id << "\"");
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  return this->values[i * this->nb_component + j];
}

template <class T, bool is_scal>
inline auto Array<T, is_scal>::operator[](Int i) -> reference {
  AKANTU_DEBUG_ASSERT(this->size_ > 0,
                      "The array \"" << this->id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < this->size_ * this->nb_component),
                      "The value at position ["
                          << i << "] is out of range in array \"" << this->id
                          << "\"");
  return this->values[i];
}

/* -------------------------------------------------------------------------- */
template <class T, bool is_scal>
inline auto Array<T, is_scal>::operator[](Int i) const -> const_reference {
  AKANTU_DEBUG_ASSERT(this->size_ > 0,
                      "The array \"" << this->id << "\" is empty");
  AKANTU_DEBUG_ASSERT((i < this->size_ * this->nb_component),
                      "The value at position ["
                          << i << "] is out of range in array \"" << this->id
                          << "\"");
  return this->values[i];
}

/* -------------------------------------------------------------------------- */
/**
 * erase an element. If the erased element is not the last of the array, the
 * last element is moved into the hole in order to maintain contiguity. This
 * may invalidate existing iterators (For instance an iterator obtained by
 * Array::end() is no longer correct) and will change the order of the
 * elements.
 * @param i index of element to erase
 */
template <class T, bool is_scal> inline void Array<T, is_scal>::erase(Int i) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT((this->size_ > 0), "The array is empty");
  AKANTU_DEBUG_ASSERT((i < this->size_), "The element at position ["
                                             << i << "] is out of range (" << i
                                             << ">=" << this->size_ << ")");

  if (i != (this->size_ - 1)) {
    for (Idx j = 0; j < this->nb_component; ++j) {
      this->values[i * this->nb_component + j] =
          this->values[(this->size_ - 1) * this->nb_component + j];
    }
  }

  this->resize(this->size_ - 1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Subtract another array entry by entry from this array in place. Both arrays
 * must
 * have the same size and nb_component. If the arrays have different shapes,
 * code compiled in debug mode will throw an expeption and optimised code
 * will behave in an unpredicted manner
 * @param other array to subtract from this
 * @return reference to modified this
 */
template <class T, bool is_scal>
Array<T, is_scal> &
Array<T, is_scal>::operator-=(const Array<T, is_scal> & vect) {
  AKANTU_DEBUG_ASSERT((this->size_ == vect.size_) &&
                          (this->nb_component == vect.nb_component),
                      "The too array don't have the same sizes");

  T * a = this->values;
  T * b = vect.data();
  for (Idx i = 0; i < this->size_ * this->nb_component; ++i) {
    *a -= *b;
    ++a;
    ++b;
  }

  return *this;
}

/* --------------------------------------------------------------------------
 */
/**
 * Add another array entry by entry to this array in
 * place. Both arrays must have the same size and
 * nb_component. If the arrays have different shapes, code
 * compiled in debug mode will throw an expeption and
 * optimised code will behave in an unpredicted manner
 * @param other array to add to this
 * @return reference to modified this
 */
template <class T, bool is_scal>
Array<T, is_scal> &
Array<T, is_scal>::operator+=(const Array<T, is_scal> & vect) {
  AKANTU_DEBUG_ASSERT((this->size_ == vect.size()) &&
                          (this->nb_component == vect.nb_component),
                      "The too array don't have the same sizes");

  T * a = this->values;
  T * b = vect.data();
  for (Idx i = 0; i < this->size_ * this->nb_component; ++i) {
    *a++ += *b++;
  }

  return *this;
}

/* --------------------------------------------------------------------------
 */
/**
 * Multiply all entries of this array by a scalar in place
 * @param alpha scalar multiplicant
 * @return reference to modified this
 */

template <class T, bool is_scal>
Array<T, is_scal> & Array<T, is_scal>::operator*=(const T & alpha) {
  T * a = this->values;
  for (Idx i = 0; i < this->size_ * this->nb_component; ++i) {
    *a++ *= alpha;
  }

  return *this;
}

/* --------------------------------------------------------------------------
 */
/**
 * Compare this array element by element to another.
 * @param other array to compare to
 * @return true it all element are equal and arrays have
 * the same shape, else false
 */
template <class T, bool is_scal>
bool Array<T, is_scal>::operator==(const Array<T, is_scal> & other) const {
  bool equal = this->nb_component == other.nb_component &&
               this->size_ == other.size_ && this->id == other.id;
  if (not equal) {
    return false;

    if (this->values == array.data()) {
      return true;
    }
    return std::equal(this->values,
                      this->values + this->size_ * this->nb_component,
                      array.data());
  }

  /* --------------------------------------------------------------------------
   */
  template <class T, bool is_scal>
  bool Array<T, is_scal>::operator!=(const Array<T, is_scal> & other) const {
    return !operator==(other);
  }

  /* --------------------------------------------------------------------------
   */
  /**
   * set all tuples of the array to a given vector or matrix
   * @param vm Matrix or Vector to fill the array with
   */
  template <class T, bool is_scal>
  template <typename C, std::enable_if_t<aka::is_tensor<C>::value> *>
  inline void Array<T, is_scal>::set(const C & vm) {
    AKANTU_DEBUG_ASSERT(
        this->nb_component == vm.size(),
        "The size of the object does not match the number of components");
    for (T * it = this->values;
         it < this->values + this->nb_component * this->size_;
         it += this->nb_component) {
      std::copy_n(vm.data(), this->nb_component, it);
    }
  }
  /* --------------------------------------------------------------------------
   */
  template <class T, bool is_scal>
  void Array<T, is_scal>::append(const Array<T> & other) {
    AKANTU_DEBUG_ASSERT(
        this->nb_component == other.nb_component,
        "Cannot append an array with a different number of component");
    Idx old_size = this->size_;
    this->resize(this->size_ + other.size());

    T * tmp = this->values + this->nb_component * old_size;
    std::copy_n(other.data(), other.size() * this->nb_component, tmp);
  }

  /* --------------------------------------------------------------------------
   */
  /* Functions Array<T, is_scal> */
  /* --------------------------------------------------------------------------
   */
  template <class T, bool is_scal>
  Array<T, is_scal>::Array(Int size, Int nb_component, const ID & id)
      : parent(size, nb_component, id) {}

  template <>
  inline Array<std::string, false>::Array(Int size, Int nb_component,
                                          const ID & id)
      : parent(size, nb_component, "", id) {}

  /* --------------------------------------------------------------------------
   */
  template <class T, bool is_scal>
  Array<T, is_scal>::Array(Int size, Int nb_component, const_reference value,
                           const ID & id)
      : parent(size, nb_component, value, id) {}

  /* --------------------------------------------------------------------------
   */
  template <class T, bool is_scal>
  Array<T, is_scal>::Array(const Array & vect, const ID & id)
      : parent(vect, id) {}

  /* --------------------------------------------------------------------------
   */
  template <class T, bool is_scal>
  Array<T, is_scal> & Array<T, is_scal>::operator=(
      const Array<T, is_scal> & other) {
    AKANTU_DEBUG_WARNING("You are copying the array "
                         << this->id << " are you sure it is on purpose");

    if (&other == this) {
      return *this;
    }

    parent::operator=(other);
    return *this;
  }

  /* --------------------------------------------------------------------------
   */
  template <class T, bool is_scal>
  Array<T, is_scal>::Array(const std::vector<T> & vect) : parent(vect) {}

  /* --------------------------------------------------------------------------
   */
  template <class T, bool is_scal> Array<T, is_scal>::~Array() = default;

  /* --------------------------------------------------------------------------
   */
  /**
   * search elem in the array, return  the position of the
   * first occurrence or -1 if not found
   *  @param elem the element to look for
   *  @return index of the first occurrence of elem or -1 if
   * elem is not present
   */
  template <class T, bool is_scal>
  Idx Array<T, is_scal>::find(const_reference elem) const {
    AKANTU_DEBUG_IN();

    auto begin = this->begin();
    auto end = this->end();
    auto it = std::find(begin, end, elem);

    AKANTU_DEBUG_OUT();
    return (it != end) ? it - begin : Idx(-1);
  }

  /* --------------------------------------------------------------------------
   */
  // template <class T, bool is_scal> UInt Array<T,
  // is_scal>::find(T elem[]) const
  // {
  //   AKANTU_DEBUG_IN();
  //   T * it = this->values;
  //   UInt i = 0;
  //   for (; i < this->size_; ++i) {
  //     if (*it == elem[0]) {
  //       T * cit = it;
  //       UInt c = 0;
  //       for (; (c < this->nb_component) && (*cit ==
  //       elem[c]); ++c, ++cit)
  //         ;
  //       if (c == this->nb_component) {
  //         AKANTU_DEBUG_OUT();
  //         return i;
  //       }
  //     }
  //     it += this->nb_component;
  //   }
  //   return UInt(-1);
  // }

  /* --------------------------------------------------------------------------
   */
  template <class T, bool is_scal>
  template <typename V, std::enable_if_t<aka::is_tensor<V>::value> *>
  inline Idx Array<T, is_scal>::find(const V & elem) {
    AKANTU_DEBUG_ASSERT(elem.size() == this->nb_component,
                        "Cannot find an element with a wrong size ("
                            << elem.size() << ") != " << this->nb_component);
    return this->find(elem.data());
  }

  /* --------------------------------------------------------------------------
   */
  /**
   * copy the content of another array. This overwrites the
   * current content.
   * @param other Array to copy into this array. It has to
   * have the same nb_component as this. If compiled in
   * debug mode, an incorrect other will result in an
   * exception being thrown. Optimised code may result in
   * unpredicted behaviour.
   * @param no_sanity_check turns off all checkes
   */
  template <class T, bool is_scal>
  void Array<T, is_scal>::copy(const Array<T, is_scal> & other,
                               bool no_sanity_check) {
    AKANTU_DEBUG_IN();

    if (not no_sanity_check and (other.nb_component != this->nb_component)) {
      AKANTU_ERROR("The two arrays do not have the same "
                   "number of components");
    }

    this->resize((other.size_ * other.nb_component) / this->nb_component);

    std::copy_n(vect.data(), this->size_ * this->nb_component, this->values);

    AKANTU_DEBUG_OUT();
  }

  /* --------------------------------------------------------------------------
   */
  template <bool is_scal> class ArrayPrintHelper {
  public:
    template <typename T>
    static void print_content(const Array<T> & vect, std::ostream & stream,
                              int indent) {
      std::string space(indent, AKANTU_INDENT);

      stream << space << " + values         : {";
      for (Idx i = 0; i < vect.size(); ++i) {
        stream << "{";
        for (Idx j = 0; j < vect.getNbComponent(); ++j) {
          stream << vect(i, j);
          if (j != vect.getNbComponent() - 1) {
            stream << ", ";
          }
        }
        stream << "}";
        if (i != vect.size() - 1) {
          stream << ", ";
        }
      }
      stream << "}" << std::endl;
    }
  };

  template <> class ArrayPrintHelper<false> {
  public:
    template <typename T>
    static void print_content(__attribute__((unused)) const Array<T> & vect,
                              __attribute__((unused)) std::ostream & stream,
                              __attribute__((unused)) int indent) {}
  };

  /* --------------------------------------------------------------------------
   */
  template <class T, bool is_scal>
  void Array<T, is_scal>::printself(std::ostream & stream, int indent) const {
    std::string space(indent, AKANTU_INDENT);

    std::streamsize prec = stream.precision();
    std::ios_base::fmtflags ff = stream.flags();

    stream.setf(std::ios_base::showbase);
    stream.precision(2);

    stream << space << "Array<" << debug::demangle(typeid(T).name()) << "> ["
           << std::endl;
    stream << space << " + id             : " << this->id << std::endl;
    stream << space << " + size           : " << this->size_ << std::endl;
    stream << space << " + nb_component   : " << this->nb_component
           << std::endl;
    stream << space << " + allocated size : " << this->getAllocatedSize()
           << std::endl;
    stream << space << " + memory size    : "
           << printMemorySize<T>(this->getMemorySize()) << std::endl;
    if (not AKANTU_DEBUG_LEVEL_IS_TEST()) {
      stream << space << " + address        : " << std::hex << this->values
             << std::dec << std::endl;
    }

    stream.precision(prec);
    stream.flags(ff);

    if (AKANTU_DEBUG_TEST(dblDump) || AKANTU_DEBUG_LEVEL_IS_TEST()) {
      ArrayPrintHelper<is_scal or std::is_enum<T>::value>::print_content(
          *this, stream, indent);
    }

    stream << space << "]" << std::endl;
  }

  /* --------------------------------------------------------------------------
   */
  /* Inline Functions ArrayBase */
  /* --------------------------------------------------------------------------
   */

  // inline bool ArrayBase::empty() { return (this->size_ ==
  // 0); }

#ifndef SWIG
  /* --------------------------------------------------------------------------
   */
  /* ArrayFilter */
  /* --------------------------------------------------------------------------
   */
  template <typename T> class ArrayFilter {
  public:
    /// const iterator
    template <class original_iterator, typename filter_iterator>
    class const_iterator {
    public:
      Idx getCurrentIndex() {
        return (*this->filter_it * this->nb_item_per_elem +
                this->sub_element_counter);
      }

      inline const_iterator() = default;
      inline const_iterator(original_iterator origin_it,
                            filter_iterator filter_it, Int nb_item_per_elem)
          : origin_it(std::move(origin_it)), filter_it(std::move(filter_it)),
            nb_item_per_elem(nb_item_per_elem), sub_element_counter(0){};

      inline bool operator!=(const_iterator & other) const {
        return !((*this) == other);
      }
      inline bool operator==(const_iterator & other) const {
        return (this->origin_it == other.origin_it &&
                this->filter_it == other.filter_it &&
                this->sub_element_counter == other.sub_element_counter);
      }

      inline bool operator!=(const const_iterator & other) const {
        return !((*this) == other);
      }
      inline bool operator==(const const_iterator & other) const {
        return (this->origin_it == other.origin_it &&
                this->filter_it == other.filter_it &&
                this->sub_element_counter == other.sub_element_counter);
      }

      inline const_iterator & operator++() {
        ++sub_element_counter;
        if (sub_element_counter == nb_item_per_elem) {
          sub_element_counter = 0;
          ++filter_it;
        }
        return *this;
      };

      inline decltype(auto) operator*() {
        return origin_it[nb_item_per_elem * (*filter_it) + sub_element_counter];
      };

    private:
      original_iterator origin_it;
      filter_iterator filter_it;
      /// the number of item per element
      Int nb_item_per_elem;
      /// counter for every sub element group
      Int sub_element_counter;
    };

    using vector_iterator = void; // iterator<Vector<T>>;
    using array_type = Array<T>;
    using const_vector_iterator =
        const_iterator<typename array_type::const_vector_iterator,
                       Array<Idx>::const_scalar_iterator>;
    using value_type = typename array_type::value_type;

  private:
    /* ------------------------------------------------------------------------
     */
    /* Constructors/Destructors */
    /* ------------------------------------------------------------------------
     */
  public:
    ArrayFilter(const Array<T> & array, const Array<Idx> & filter,
                Int nb_item_per_elem)
        : array(array), filter(filter), nb_item_per_elem(nb_item_per_elem){};

    decltype(auto) begin_reinterpret(Int n, Int new_size) const {
      Int new_nb_item_per_elem = this->nb_item_per_elem;
      if (new_size != 0 && n != 0)
        new_nb_item_per_elem = this->array.getNbComponent() *
                               this->filter.size() * this->nb_item_per_elem /
                               (n * new_size);

      return const_vector_iterator(make_view(this->array, n).begin(),
                                   this->filter.begin(), new_nb_item_per_elem);
    };

    decltype(auto) end_reinterpret(Int n, Int new_size) const {
      Int new_nb_item_per_elem = this->nb_item_per_elem;
      if (new_size != 0 && n != 0)
        new_nb_item_per_elem = this->array.getNbComponent() *
                               this->filter.size() * this->nb_item_per_elem /
                               (n * new_size);

      return const_vector_iterator(make_view(this->array, n).begin(),
                                   this->filter.end(), new_nb_item_per_elem);
    };

    // vector_iterator begin_reinterpret(Int, Int) { throw; };
    // vector_iterator end_reinterpret(Int, Int) { throw; };

    /// return the size of the filtered array which is the filter size
    Int size() const { return this->filter.size() * this->nb_item_per_elem; };
    /// the number of components of the filtered array
    Int getNbComponent() const { return this->array.getNbComponent(); };

    /* ------------------------------------------------------------------------
     */
    /* Class Members */
    /* ------------------------------------------------------------------------
     */
  private:
    /// reference to array of data
    const Array<T> & array;
    /// reference to the filter used to select elements
    const Array<Idx> & filter;
    /// the number of item per element
    Int nb_item_per_elem;
  };

  /* --------------------------------------------------------------------------
   */
  /* Begin/End functions implementation */
  /* --------------------------------------------------------------------------
   */
  namespace detail {
    template <class Tuple, size_t... Is>
    constexpr auto take_front_impl(Tuple && t,
                                   std::index_sequence<Is...> /*idxs*/) {
      return std::make_tuple(std::get<Is>(std::forward<Tuple>(t))...);
    }

    template <size_t N, class Tuple> constexpr auto take_front(Tuple && t) {
      return take_front_impl(std::forward<Tuple>(t),
                             std::make_index_sequence<N>{});
    }

    template <typename... T> std::string to_string_all(T &&... t) {
      if (sizeof...(T) == 0) {
        return "";
      }

      std::stringstream ss;
      bool noComma = true;
      ss << "(";
      (void)std::initializer_list<bool>{
          (ss << (noComma ? "" : ", ") << t, noComma = false)...};
      ss << ")";
      return ss.str();
    }

    template <std::size_t N> struct InstantiationHelper {
      template <typename type, typename T, typename... Ns>
      static auto instantiate(T && data, Ns... ns) {
        return std::make_unique<type>(data, ns...);
      }
    };

    template <> struct InstantiationHelper<0> {
      template <typename type, typename T> static auto instantiate(T && data) {
        return data;
      }
    };

    template <typename Arr, typename T, typename... Ns>
    decltype(auto) get_iterator(Arr && array, T * data, Ns &&... ns) {
      const auto is_const_arr =
          std::is_const<std::remove_reference_t<Arr>>::value;
      using type = ViewIteratorHelper_t<sizeof...(Ns) - 1, T>;
      using iterator =
          std::conditional_t<is_const_arr, const_view_iterator<type>,
                             view_iterator<type>>;
      static_assert(sizeof...(Ns), "You should provide a least one size");

      if (array.getNbComponent() * array.size() !=
          Int(product_all(std::forward<Ns>(ns)...))) {
        AKANTU_CUSTOM_EXCEPTION_INFO(
            debug::ArrayException(),
            "The iterator on "
                << debug::demangle(typeid(Arr).name())
                << to_string_all(array.size(), array.getNbComponent())
                << "is not compatible with the type "
                << debug::demangle(typeid(type).name())
                << to_string_all(ns...));
      }

      auto && it =
          aka::apply([&](auto... n) { return iterator(data, n...); },
                     take_front<sizeof...(Ns) - 1>(std::make_tuple(ns...)));

      return it;
    }
  } // namespace detail

  /* --------------------------------------------------------------------------
   */
  /* --------------------------------------------------------------------------
   */
  template <class T, bool is_scal>
  template <typename... Ns>
  inline auto Array<T, is_scal>::begin(Ns && ... ns) {
    return detail::get_iterator(*this, this->values, std::forward<Ns>(ns)...,
                                this->size_);
  }

  template <class T, bool is_scal>
  template <typename... Ns>
  inline auto Array<T, is_scal>::end(Ns && ... ns) {
    return detail::get_iterator(*this,
                                this->values + this->nb_component * this->size_,
                                std::forward<Ns>(ns)..., this->size_);
  }

  template <class T, bool is_scal>
  template <typename... Ns>
  inline auto Array<T, is_scal>::begin(Ns && ... ns) const {
    return detail::get_iterator(*this, this->values, std::forward<Ns>(ns)...,
                                this->size_);
  }

  template <class T, bool is_scal>
  template <typename... Ns>
  inline auto Array<T, is_scal>::end(Ns && ... ns) const {
    return detail::get_iterator(*this,
                                this->values + this->nb_component * this->size_,
                                std::forward<Ns>(ns)..., this->size_);
  }

  template <class T, bool is_scal>
  template <typename... Ns>
  inline auto Array<T, is_scal>::begin_reinterpret(Ns && ... ns) {
    return detail::get_iterator(*this, this->values, std::forward<Ns>(ns)...);
  }

  template <class T, bool is_scal>
  template <typename... Ns>
  inline auto Array<T, is_scal>::end_reinterpret(Ns && ... ns) {
    return detail::get_iterator(
        *this, this->values + detail::product_all(std::forward<Ns>(ns)...),
        std::forward<Ns>(ns)...);
  }

  template <class T, bool is_scal>
  template <typename... Ns>
  inline auto Array<T, is_scal>::begin_reinterpret(Ns && ... ns) const {
    return detail::get_iterator(*this, this->values, std::forward<Ns>(ns)...);
  }

  template <class T, bool is_scal>
  template <typename... Ns>
  inline auto Array<T, is_scal>::end_reinterpret(Ns && ... ns) const {
    return detail::get_iterator(
        *this, this->values + detail::product_all(std::forward<Ns>(ns)...),
        std::forward<Ns>(ns)...);
  }

  /* --------------------------------------------------------------------------
   */
  /* Views */
  /* --------------------------------------------------------------------------
   */
  namespace detail {
    template <typename Array, typename... Ns> class ArrayView {
      using tuple = std::tuple<Ns...>;

    public:
      using size_type = Idx;

      ~ArrayView() = default;
      constexpr ArrayView(Array && array, Ns... ns) noexcept
          : array(array), sizes(std::move(ns)...) {}

      constexpr ArrayView(const ArrayView & array_view) = default;
      constexpr ArrayView(ArrayView && array_view) noexcept = default;

      constexpr ArrayView & operator=(const ArrayView & array_view) = default;
      constexpr ArrayView &
      operator=(ArrayView && array_view) noexcept = default;

      auto begin() {
        return aka::apply(
            [&](auto &&... ns) {
              return detail::get_iterator(array.get(), array.get().data(),
                                          std::forward<decltype(ns)>(ns)...);
            },
            sizes);
      }

      auto begin() const {
        return aka::apply(
            [&](auto &&... ns) {
              return detail::get_iterator(array.get(), array.get().data(),
                                          std::forward<decltype(ns)>(ns)...);
            },
            sizes);
      }

      auto end() {
        return aka::apply(
            [&](auto &&... ns) {
              return detail::get_iterator(
                  array.get(),
                  array.get().data() +
                      detail::product_all(std::forward<decltype(ns)>(ns)...),
                  std::forward<decltype(ns)>(ns)...);
            },
            sizes);
      }

      auto end() const {
        return aka::apply(
            [&](auto &&... ns) {
              return detail::get_iterator(
                  array.get(),
                  array.get().data() +
                      detail::product_all(std::forward<decltype(ns)>(ns)...),
                  std::forward<decltype(ns)>(ns)...);
            },
            sizes);
      }

      auto cbegin() const { return this->begin(); }
      auto cend() const { return this->end(); }

      constexpr auto size() const {
        return std::get<std::tuple_size<tuple>::value - 1>(sizes);
      }

      constexpr auto dims() const { return std::tuple_size<tuple>::value - 1; }

    private:
      std::reference_wrapper<std::remove_reference_t<Array>> array;
      tuple sizes;
    };

    /* ------------------------------------------------------------------------
     */
    template <typename T, typename... Ns>
    class ArrayView<const ArrayFilter<T> &, Ns...> {
      using tuple = std::tuple<Ns...>;

    public:
      constexpr ArrayView(const ArrayFilter<T> & array, Ns... ns)
          : array(array), sizes(std::move(ns)...) {}

      constexpr ArrayView(const ArrayView & array_view) = default;
      constexpr ArrayView(ArrayView && array_view) = default;

      constexpr ArrayView & operator=(const ArrayView & array_view) = default;
      constexpr ArrayView & operator=(ArrayView && array_view) = default;

      auto begin() const {
        return aka::apply(
            [&](auto &&... ns) {
              return array.get().begin_reinterpret(
                  std::forward<decltype(ns)>(ns)...);
            },
            sizes);
      }

      auto end() const {
        return aka::apply(
            [&](auto &&... ns) {
              return array.get().end_reinterpret(
                  std::forward<decltype(ns)>(ns)...);
            },
            sizes);
      }

      auto cbegin() const { return this->begin(); }
      auto cend() const { return this->end(); }

      constexpr auto size() const {
        return std::get<std::tuple_size<tuple>::value - 1>(sizes);
      }

      constexpr auto dims() const { return std::tuple_size<tuple>::value - 1; }

    private:
      std::reference_wrapper<const ArrayFilter<T>> array;
      tuple sizes;
    };
  } // namespace detail

  /* --------------------------------------------------------------------------
   */
  template <typename Array, typename... Ns>
  decltype(auto) make_view(Array && array, const Ns... ns) {
    static_assert(
        aka::conjunction<std::is_integral<std::decay_t<Ns>>...>::value,
        "Ns should be integral types");
    AKANTU_DEBUG_ASSERT((detail::product_all(ns...) != 0),
                        "You must specify non zero dimensions");
    auto size = std::forward<decltype(array)>(array).size() *
                std::forward<decltype(array)>(array).getNbComponent() /
                detail::product_all(ns...);

    return detail::ArrayView<Array, std::common_type_t<size_t, Ns>...,
                             std::common_type_t<size_t, decltype(size)>>(
        std::forward<Array>(array), std::move(ns)..., size);
  }

  /* --------------------------------------------------------------------------
   */
  template <class T, bool is_scal>
  template <typename R>
  inline auto Array<T, is_scal>::erase(const view_iterator<R> & it) {
    T * curr = it.data();
    Idx pos = (curr - this->values) / this->nb_component;
    erase(pos);
    view_iterator<R> rit = it;
    return --rit;
  }

} // namespace akantu

//#endif /* __AKANTU_AKA_ARRAY_TMPL_HH__ */
