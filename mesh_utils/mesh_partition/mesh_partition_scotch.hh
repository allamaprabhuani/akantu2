/**
 * @file   mesh_partition_scotch.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Aug 13 10:00:06 2010
 *
 * @brief  mesh partitioning based on libScotch
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MESH_PARTITION_SCOTCH_HH__
#define __AKANTU_MESH_PARTITION_SCOTCH_HH__

/* -------------------------------------------------------------------------- */
#include "mesh_partition.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class MeshPartitionScotch : public MeshPartition {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MeshPartitionScotch();

  virtual ~MeshPartitionScotch();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:



  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
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

#include "mesh_partition_scotch_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MeshPartitionScotch & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_MESH_PARTITION_SCOTCH_HH__ */
