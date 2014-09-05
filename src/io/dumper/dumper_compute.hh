#ifndef __AKANTU_DUMPER_COMPUTE_HH__
#define __AKANTU_DUMPER_COMPUTE_HH__
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "dumper_iohelper.hh"
#include "dumper_type_traits.hh"
#include "dumper_field.hh"
#include "io_helper.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__

class ComputeFunctorInterface {

public:

  virtual ~ComputeFunctorInterface(){};

  virtual UInt getDim() = 0;
  virtual UInt getNbComponent(UInt old_nb_comp) = 0;

};

/* -------------------------------------------------------------------------- */

template <typename return_type>
class ComputeFunctorOutput : public ComputeFunctorInterface {

public:

  ComputeFunctorOutput(){};
  virtual ~ComputeFunctorOutput(){};

};


/* -------------------------------------------------------------------------- */

template <typename input_type,typename return_type>
class ComputeFunctor : public ComputeFunctorOutput<return_type> {

public:

  ComputeFunctor(){};
  virtual ~ComputeFunctor(){};

  virtual return_type func(const input_type & d, Element global_index) = 0;
};

/* -------------------------------------------------------------------------- */


template <typename SubFieldCompute, typename _return_type>
class FieldCompute : public Field {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */  
  
public:
 
  typedef typename SubFieldCompute::iterator sub_iterator;
  typedef typename SubFieldCompute::types sub_types;
  typedef typename sub_types::return_type sub_return_type;
  typedef _return_type return_type;
  typedef typename sub_types::data_type data_type;

  typedef TypeTraits<data_type, return_type, ElementTypeMapArray<data_type> > types;


  class iterator {
  public:
    iterator(const sub_iterator & it, ComputeFunctor<sub_return_type,return_type> & func) 
      : it(it), func(func) {}
    
    bool operator!=(const iterator & it)const { return it.it != this->it; }
    iterator operator++() { ++this->it; return *this; }
    

    UInt currentGlobalIndex(){
      return this->it.currentGlobalIndex();
    }

    return_type operator*() {
      return func.func(*it,it.getCurrentElement());
    }

    Element getCurrentElement(){
      return this->it.getCurrentElement();
    }

    
  protected:
    sub_iterator it;
    ComputeFunctor<sub_return_type,return_type> & func;
  };


  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:

  FieldCompute(SubFieldCompute & cont, ComputeFunctorInterface & func)
    :sub_field(cont),func(dynamic_cast<ComputeFunctor<sub_return_type,return_type> &>(func)){
    this->checkHomogeneity();
  };

  virtual void registerToDumper(const std::string & id,
				iohelper::Dumper & dumper) {
    dumper.addElemDataField(id, *this);
  }


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

  iterator begin() { return iterator(sub_field.begin(), func); }
  iterator end  () { return iterator(sub_field.end(),   func); }


  UInt getDim() {
    return func.getDim();
  }


  UInt size() { 
    throw;
    // return Functor::size();
    return 0;
  }

  virtual void checkHomogeneity(){this->homogeneous = true;};

  iohelper::DataType getDataType() { return iohelper::getDataType<data_type>(); }


  /// get the number of components of the hosted field
  virtual ElementTypeMap<UInt> getNbComponents(UInt dim = _all_dimensions,
					      GhostType ghost_type = _not_ghost,
					      ElementKind kind = _ek_not_defined){

    ElementTypeMap<UInt> nb_components;
    const ElementTypeMap<UInt> & old_nb_components = 
      this->sub_field.getNbComponents(dim,ghost_type,kind);
    

    ElementTypeMap<UInt>::type_iterator tit = old_nb_components.firstType(dim,ghost_type,kind);
    ElementTypeMap<UInt>::type_iterator end = old_nb_components.lastType(dim,ghost_type,kind);

    while (tit != end){
      UInt nb_comp = old_nb_components(*tit,ghost_type);
      nb_components(*tit,ghost_type) = func.getNbComponent(nb_comp);
      ++tit;
    }
  return nb_components;
};

  /// for connection to a FieldCompute
  inline virtual Field * connect(FieldComputeProxy & proxy);

  /// for connection to a FieldCompute
  virtual ComputeFunctorInterface * connect(HomogenizerProxy & proxy);


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

public:

  SubFieldCompute & sub_field;
  ComputeFunctor<sub_return_type,return_type> & func;
  
};


/* -------------------------------------------------------------------------- */

class FieldComputeProxy {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:

  FieldComputeProxy(ComputeFunctorInterface & func):func(func){};

public:

  inline static Field * createFieldCompute(Field * field, 
					   ComputeFunctorInterface & func){
    
    FieldComputeProxy compute_proxy(func);
    return field->connect(compute_proxy);
  }

  template <typename T>
  Field * connectToField(T * ptr){
    if (dynamic_cast<ComputeFunctorOutput<Vector<Real> > *>(&func)){
      return this->connectToFunctor<Vector<Real> >(ptr);
    }
    else if (dynamic_cast<ComputeFunctorOutput<Vector<UInt> > *>(&func)){
      return this->connectToFunctor<Vector<UInt> >(ptr);
    }

    else if (dynamic_cast<ComputeFunctorOutput<Matrix<UInt> > *>(&func)){
      return this->connectToFunctor<Matrix<UInt> >(ptr);
    }

    else if (dynamic_cast<ComputeFunctorOutput<Matrix<Real> > *>(&func)){
      return this->connectToFunctor<Matrix<Real> >(ptr);
    }


    else throw;
  }

  template <typename output, typename T>
  Field * connectToFunctor(T * ptr){
    return new FieldCompute<T,output>(*ptr,func);
  }

  template <typename output, typename SubFieldCompute, 
	    typename return_type1, typename return_type2>
  Field * connectToFunctor(FieldCompute<FieldCompute<SubFieldCompute,return_type1>,
			   return_type2> * ptr){
    throw;    //    return new FieldCompute<T,output>(*ptr,func);
    return NULL;
  }


  template <typename output, typename SubFieldCompute, 
   	    typename return_type1,typename return_type2,
   	    typename return_type3,typename return_type4>
  Field * connectToFunctor(FieldCompute<
    			     FieldCompute<
   			       FieldCompute<
   			         FieldCompute<SubFieldCompute,return_type1>,
   			           return_type2>,
   			         return_type3>,
   			       return_type4> * ptr){
     throw;    //    return new FieldCompute<T,output>(*ptr,func);
     return NULL;
   }
            

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

  ComputeFunctorInterface & func;
 
};


/* -------------------------------------------------------------------------- */

/// for connection to a FieldCompute
template <typename SubFieldCompute, typename return_type>
inline Field * FieldCompute<SubFieldCompute,return_type>
::connect(FieldComputeProxy & proxy){

  return proxy.connectToField(this);
}

/* -------------------------------------------------------------------------- */






__END_AKANTU_DUMPER__
__END_AKANTU__


#endif /* __AKANTU_DUMPER_COMPUTE_HH__ */
