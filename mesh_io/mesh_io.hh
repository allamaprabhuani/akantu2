/**
 * @file   mesh_io.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jun 18 10:27:42 2010
 *
 * @brief  Interface of a mesh io class, reader and writer
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MESH_IO_HH__
#define __AKANTU_MESH_IO_HH__

/* -------------------------------------------------------------------------- */
#include "common.hh"
#include "mesh.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class MeshIO {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MeshIO();

  virtual ~MeshIO();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {};

  /// read a mesh from the file
  virtual void read(const std::string & filename, const Mesh & mesh) = 0;

  /// write a mesh to a file
  virtual void write(const std::string & filename, const Mesh & mesh) = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  bool canReadSurface;

  bool canReadExtendedData;

  //  std::string filename;

  //  Mesh & mesh;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MeshIO & _this)
{
  _this.printself(stream);

  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_MESH_IO_HH__ */
