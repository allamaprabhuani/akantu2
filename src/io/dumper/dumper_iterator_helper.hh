#ifndef __AKANTU_DUMPER_ITERATOR_HELPER_HH__
#define __AKANTU_DUMPER_ITERATOR_HELPER_HH__
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__

template<class T, class R>
class iterator_helper {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */  
public:
  

  typedef typename Array<T>::template const_iterator< R > internal_iterator;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  
public:

  static internal_iterator begin(const Array<T> & vect, UInt m, UInt n, UInt size) {
    return vect.begin_reinterpret(n*m, size);
  }

  static internal_iterator end(const Array<T> & vect, UInt m, UInt n, UInt size) {
    return vect.end_reinterpret(n*m, size);
  }
};

/* -------------------------------------------------------------------------- */

template<class T>
class iterator_helper<T, Matrix<T> > {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */  

public:

  typedef typename Array<T>::template const_iterator< Matrix<T> > internal_iterator;

public:

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  

  static internal_iterator begin(const Array<T> & vect, UInt m, UInt n, UInt size) {
    return vect.begin_reinterpret(m, n, size);
  }

  static internal_iterator end(const Array<T> & vect, UInt m, UInt n, UInt size) {
    return vect.end_reinterpret(m, n, size);
  }
};


/* -------------------------------------------------------------------------- */

__END_AKANTU_DUMPER__
__END_AKANTU__


#endif /* __AKANTU_DUMPER_ITERATOR_HELPER_HH__ */
