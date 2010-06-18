/**
 * @file   mesh_io_msh.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jun 18 11:30:59 2010
 *
 * @brief  Read/Write for MSH files
 *
 * @section LICENSE
 *
 * <insert lisence here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __MYFEM_MESH_IO_MSH_HH__
#define __MYFEM_MESH_IO_MSH_HH__

/* -------------------------------------------------------------------------- */
#include "mesh_io.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_MYFEM__

class MeshIOMSH {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MeshIOMSH();
  virtual ~MeshIOMSH();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// read a mesh from the file
  static virtual void read(const std::string & filename, const Mesh & mesh);

  /// write a mesh to a file
  static virtual void write(const std::string & filename, const Mesh & mesh);

  /* ------------------------------------------------------------------------ */
  /* Accesors                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MeshIOMSH & _this)
{
  _this.printself(stream);
}


__END_MYFEM__

#endif /* __MYFEM_MESH_IO_MSH_HH__ */
