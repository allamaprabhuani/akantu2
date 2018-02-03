/**
 * @file   ntn_fricreg_rubin_ampuero.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  regularisation that regularizes the contact pressure
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __AST_NTN_FRICREG_RUBIN_AMPUERO_HH__
#define __AST_NTN_FRICREG_RUBIN_AMPUERO_HH__

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_fricreg_no_regularisation.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
class NTNFricRegRubinAmpuero : public NTNFricRegNoRegularisation {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  NTNFricRegRubinAmpuero(NTNBaseContact * contact,
			 const FrictionID & id = "rubin_ampuero",
			 const MemoryID & memory_id = 0);
  virtual ~NTNFricRegRubinAmpuero() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void registerSynchronizedArray(SynchronizedArrayBase & array);
  virtual void dumpRestart(const std::string & file_name) const;
  virtual void readRestart(const std::string & file_name);

  virtual void setToSteadyState();

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  virtual void addDumpFieldToDumper(const std::string & dumper_name,
				    const std::string & field_id);
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

protected:
  /// get the contact pressure (the norm: scalar value)
  virtual const SynchronizedArray<Real> & internalGetContactPressure();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  SynchronizedArray<Real> t_star;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ntn_fricreg_rubin_ampuero_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, 
				  const NTNFricRegRubinAmpuero & _this)
{
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* __AST_NTN_FRICREG_RUBIN_AMPUERO_HH__ */
