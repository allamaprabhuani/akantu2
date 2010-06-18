/**
 * @file   array_info.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Jun 14 17:58:46 2010
 *
 * @brief
 *
 * @section LICENSE
 *
 * <insert lisence here>
 *
 */


/* -------------------------------------------------------------------------- */
#ifndef __MYFEM_ARRAY_INFO_HH__
#define __MYFEM_ARRAY_INFO_HH__

/* -------------------------------------------------------------------------- */
#include "common.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_MYFEM__

/**
 * @class ArrayInfo
 * @brief descriptive information for arrays
 */
class ArrayInfo {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /// empty constructor
  ArrayInfo() {};

  /// complete constructor
  ArrayInfo(std::string name,
	    unsigned int nb_tuples,
	    unsigned int nb_component,
	    TypeCode type,
	    void * adresse) {
    this->name = name;
    this->nb_reserved_tuples = nb_tuples;
    this->nb_tuples = nb_tuples;
    this->nb_component = nb_component;
    this->adresse = adresse;
    this->type = type;
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  MYFEM_SET_MACRO(NbReservedTuples, nb_reserved_tuples, unsigned int);
  MYFEM_GET_MACRO(NbReservedTuples, nb_reserved_tuples, unsigned int);

  MYFEM_SET_MACRO(NbTuples, nb_tuples, unsigned int);
  MYFEM_GET_MACRO(NbTuples, nb_tuples, unsigned int);

  MYFEM_GET_MACRO(Name, name, std::string);

  MYFEM_GET_MACRO(NbComponent, nb_component, unsigned int);

  MYFEM_SET_MACRO(Adresse, adresse, void *);
  MYFEM_GET_MACRO(Adresse, adresse, void *);

  MYFEM_GET_MACRO(Type, type, TypeCode);

  /// get the size of memory allocated
  unsigned int getMemorySize() {
    return nb_reserved_tuples * nb_component * MYFEM_SIZEOF(type);
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// the name of the array
  std::string name;

  /// the number of tuples allocated
  unsigned int nb_reserved_tuples;

  /// the number of tuples used
  unsigned int nb_tuples;

  /// number of components
  unsigned int nb_component;

  /// memory adresse of the array
  void *adresse;

  /// type of the array, see common.hpp to have the list of types
  TypeCode type;
};

__END_MYFEM__

#endif /* __MYFEM_ARRAY_INFO_HH__ */
