#ifndef __AKANTU_DUMPER_FIELD_HH__
#define __AKANTU_DUMPER_FIELD_HH__
/* -------------------------------------------------------------------------- */
#include "dumper_iohelper.hh"
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__
/* -------------------------------------------------------------------------- */
class FieldComputeProxy;
class FieldComputeBaseInterface;
class ComputeFunctorInterface;
class HomogenizerProxy;
/* -------------------------------------------------------------------------- */


/// Field interface
class Field {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:

  Field(): homogeneous(false) {}
  virtual ~Field() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  /// register this to the provided dumper
  virtual void registerToDumper(const std::string & id, 
				iohelper::Dumper & dumper) = 0;

  
  /// set the number of data per item (used for elements fields at the moment)
  virtual void setNbData(UInt nb_data){AKANTU_DEBUG_TO_IMPLEMENT();};

  /// set the number of data per elem (used for elements fields at the moment)
  virtual void setNbDataPerElem(const ElementTypeMap<UInt> & nb_data){AKANTU_DEBUG_TO_IMPLEMENT();};

  /// set the number of data per elem (used for elements fields at the moment)
  virtual void setNbDataPerElem(UInt nb_data){AKANTU_DEBUG_TO_IMPLEMENT();};

  /// get the number of components of the hosted field
  virtual ElementTypeMap<UInt> getNbComponents(){throw;};

  /// for connection to a FieldCompute
  inline virtual Field * connect(FieldComputeProxy & proxy){throw;};

  /// for connection to a FieldCompute
  inline virtual ComputeFunctorInterface * connect(HomogenizerProxy & proxy){throw;};


  /// check if the same quantity of data for all element types
  virtual void checkHomogeneity() = 0;

  /// return the dumper name
  std::string getGroupName(){return group_name;};

  /// return the id of the field
  std::string getID(){return field_id;};

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
  
  /// return the flag to know if the field is homogeneous/contiguous
  virtual bool isHomogeneous() { return homogeneous; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  
protected:
  /// the flag to know if it is homogeneous
  bool homogeneous;

  /// the name of the group it was associated to
  std::string group_name;

  /// the name of the dumper it was associated to
  std::string field_id;
};

/* -------------------------------------------------------------------------- */






__END_AKANTU_DUMPER__
__END_AKANTU__


#endif /* __AKANTU_DUMPER_FIELD_HH__ */
