/**
 * @file   mesh_io_msh.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jun 18 11:30:59 2010
 *
 * @brief  Read/Write for MSH files
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MESH_IO_MSH_HH__
#define __AKANTU_MESH_IO_MSH_HH__

/* -------------------------------------------------------------------------- */
#include "mesh_io.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__


class MeshIOMSH : public MeshIO {
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
  virtual void printself(std::ostream & stream, int indent = 0) const {};

  /// read a mesh from the file
  virtual void read(const std::string & filename, const Mesh & mesh);

  /// write a mesh to a file
  virtual void write(const std::string & filename, const Mesh & mesh);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// MSH element types
  enum MSHElementType {
    _msh_not_defined   = 0,
    _msh_line_1        = 1,   // 2-node line.
    _msh_triangle_1    = 2,   // 3-node triangle.
    _msh_quadrangle_1  = 3,   // 4-node quadrangle.
    _msh_tetrahedron_1 = 4,   // 4-node tetrahedron.
    _msh_hexaedron_1   = 5,   // 8-node hexahedron.
    _msh_prism_1       = 6,   // 6-node prism.
    _msh_pyramid_1     = 7,   // 5-node pyramid.
    _msh_line_2        = 8,   // 3-node second order line
    _msh_triangle_2    = 9,   // 6-node second order triangle
    _msh_quadrangle_2  = 10,  // 9-node second order quadrangle
    _msh_tetrahedron_2 = 11,  // 10-node second order tetrahedron
    _msh_hexaedron_2   = 12,  // 27-node second order hexahedron
    _msh_prism_2       = 13,  // 18-node second order prism
    _msh_pyramid_2     = 14,  // 14-node second order pyramid
    _msh_point         = 15   // 1-node point.
  };


  /// order in witch element as to be read
  static unsigned int _read_order[_max_element_type][MAX_NUMBER_OF_NODE_PER_ELEMENT];

  /// number of nodes per msh element
  static unsigned int _msh_nodes_per_elem[16]; // 16 = number of recognized
                                               // msh element types +1 (for 0)

  /// correspondance between msh element types and akantu element types
  static ElementType _msh_to_akantu_element_types[16];

  /// correspondance between akantu element types and msh element types
  static MSHElementType _akantu_to_msh_element_types[_max_element_type];
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MeshIOMSH & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_MESH_IO_MSH_HH__ */
