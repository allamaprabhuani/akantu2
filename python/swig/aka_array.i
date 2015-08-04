%{
#define SWIG_FILE_WITH_INIT
#include "aka_array.hh"
%}

%include "typemaps.i"

namespace akantu {
  %ignore Array::operator=;
  %ignore Array::operator[];
  %ignore Array::operator();
  %ignore Array::set;
  %ignore Array::begin;
  %ignore Array::end;
  %ignore Array::begin_reinterpret;
  %ignore Array::end_reinterpret;
};

%include "aka_array.hh"

namespace akantu {
  %template(RArray) Array<akantu::Real, true>;
  %template(UArray) Array<akantu::UInt, true>;
  %template(BArray) Array<bool, true>;
}



%include "numpy.i"
%init %{
  import_array();
%}

%inline %{
  namespace akantu{
    template <typename T>
      class ArrayForPython : public Array<T>{

    public:
    ArrayForPython(T * wrapped_memory,
		   UInt size = 0,
		   UInt nb_component = 1,
		   const ID & id = "")
      : Array<T>(0,nb_component,id){
	this->values = wrapped_memory;
	this->size = size;
      };

      ~ArrayForPython(){
	this->values = NULL;
      };

      void resize(UInt new_size){
	AKANTU_DEBUG_ASSERT(this->size == new_size,"cannot resize a temporary vector");
      }
    
    };
  }
 template <typename T> int getPythonDataTypeCode(){ AKANTU_EXCEPTION("undefined type");}
 template <> int getPythonDataTypeCode<bool>(){
   int data_typecode = NPY_NOTYPE;
   size_t s = sizeof(bool);
   switch(s) {
   case 1: data_typecode = NPY_BOOL; break;
   case 2: data_typecode = NPY_UINT16; break;
   case 4: data_typecode = NPY_UINT32; break;
   case 8: data_typecode = NPY_UINT64; break;
   }
   return data_typecode;
 }

 template <> int getPythonDataTypeCode<double>(){return NPY_DOUBLE;}
 template <> int getPythonDataTypeCode<long double>(){return NPY_LONGDOUBLE;}
 template <> int getPythonDataTypeCode<float>(){return NPY_FLOAT;}
 template <> int getPythonDataTypeCode<unsigned long>(){
   int data_typecode = NPY_NOTYPE;
   size_t s = sizeof(unsigned long);
   switch(s) {
   case 2: data_typecode  = NPY_UINT16;  break;
   case 4: data_typecode  = NPY_UINT32; break;
   case 8: data_typecode  = NPY_UINT64; break;
   }
   return data_typecode;
 }
 template <> int getPythonDataTypeCode<akantu::UInt>(){
   int data_typecode = NPY_NOTYPE;
   size_t s = sizeof(akantu::UInt);
   switch(s) {
   case 2: data_typecode  = NPY_UINT16;  break;
   case 4: data_typecode  = NPY_UINT32; break;
   case 8: data_typecode  = NPY_UINT64; break;
   }
   return data_typecode;
 }
 template <> int getPythonDataTypeCode<int>(){
   int data_typecode = NPY_NOTYPE;
   size_t s = sizeof(int);
   switch(s) {
   case 2: data_typecode = NPY_INT16; break;
   case 4: data_typecode = NPY_INT32; break;
   case 8: data_typecode = NPY_INT64; break;
   }
   return data_typecode;
 }

 int getSizeOfPythonType(int type_num){
   switch (type_num){
   case NPY_INT16 : return 2;break;
   case NPY_UINT16: return 2;break;
   case NPY_INT32 : return 4;break;
   case NPY_UINT32: return 4;break;
   case NPY_INT64 : return 8;break;
   case NPY_UINT64: return 8;break;
   case NPY_FLOAT:  return sizeof(float);break;
   case NPY_DOUBLE:  return sizeof(double);break;
   case NPY_LONGDOUBLE:  return sizeof(long double);break;
   }
   return 0;
 }

 std::string getPythonTypeName(int type_num){
   switch (type_num){
   case NPY_INT16 : return "NPY_INT16" ;break;
   case NPY_UINT16: return "NPY_UINT16";break;
   case NPY_INT32 : return "NPY_INT32" ;break;
   case NPY_UINT32: return "NPY_UINT32";break;
   case NPY_INT64 : return "NPY_INT64" ;break;
   case NPY_UINT64: return "NPY_UINT64";break;
   case NPY_FLOAT:  return "NPY_FLOAT" ;break;
   case NPY_DOUBLE: return "NPY_DOUBLE";break;
   case NPY_LONGDOUBLE: return "NPY_LONGDOUBLE";break;
   }
   return 0;
 }

 template <typename T> void checkDataType(int type_num){
   AKANTU_DEBUG_ASSERT(type_num == getPythonDataTypeCode<T>(),
		       "incompatible types between numpy and input function: " <<
		       type_num << " != " << getPythonDataTypeCode<T>() << std::endl <<
		       getSizeOfPythonType(type_num) << " != " << sizeof(T) <<
		       std::endl <<
		       "The numpy array is of type " << getPythonTypeName(type_num));
 }

 %}


%define %akantu_array_typemaps(DATA_TYPE)

%typemap(out, fragment="NumPy_Fragments") (akantu::Array< DATA_TYPE > &)
{

  int data_typecode = getPythonDataTypeCode< DATA_TYPE >();
  npy_intp dims[2] = {$1->getSize(), $1->getNbComponent()};
  PyObject* obj = PyArray_SimpleNewFromData(2, dims, data_typecode, $1->storage());
  PyArrayObject* array = (PyArrayObject*) obj;
  if (!array) SWIG_fail;
  $result = SWIG_Python_AppendOutput($result, obj);
}

%typemap(in) akantu::Array< DATA_TYPE > &
{
  if (!PyArray_Check($input)) {
    AKANTU_EXCEPTION("incompatible input which is not a numpy");
  }
  else {
    PyArray_Descr * numpy_type = (PyArray_Descr*)PyArray_DESCR((PyArrayObject*)$input);
    int type_num = numpy_type->type_num;
    checkDataType< DATA_TYPE >(type_num);
    UInt _n = PyArray_NDIM((PyArrayObject*)$input);
    if (_n != 2) AKANTU_EXCEPTION("incompatible numpy dimension " << _n);
    npy_intp * ndims = PyArray_DIMS((PyArrayObject*)$input);
    akantu::UInt sz = ndims[0];
    akantu::UInt nb_components = ndims[1];
    PyArrayIterObject *iter = (PyArrayIterObject *)PyArray_IterNew($input);
    if (iter == NULL) AKANTU_EXCEPTION("Python internal error");
    $1 = new akantu::ArrayForPython< DATA_TYPE >((DATA_TYPE*)(iter->dataptr),sz,nb_components,"tmp_array_for_python");
  }
}
%enddef


%akantu_array_typemaps(double      )
%akantu_array_typemaps(float       )
%akantu_array_typemaps(unsigned int)
%akantu_array_typemaps(unsigned long)
%akantu_array_typemaps(int         )
%akantu_array_typemaps(bool        )
