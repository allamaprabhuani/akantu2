/**
 * @file   dof_synchronizer.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon May 23 16:50:47 2011
 *
 * @brief  Synchronize Vector of DOFs
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_vector.hh"
#include "static_communicator.hh"

/* -------------------------------------------------------------------------- */



#ifndef __AKANTU_DOF_SYNCHRONIZER_HH__
#define __AKANTU_DOF_SYNCHRONIZER_HH__

__BEGIN_AKANTU__

class Mesh;

template<typename T>
class AddOperation {
public:
  inline T operator()(T & b, T & a) { return a + b; };
};


class DOFSynchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  DOFSynchronizer(const Mesh & mesh, UInt nb_degree_of_freedom);
  virtual ~DOFSynchronizer();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// init the scheme for scatter and gather operation, need extra memory
  void initScatterGatherCommunicationScheme();

  /// initialize the equation number with local ids
  void initLocalDOFEquationNumbers();

  /// initialize the equation number with local ids
  void initGlobalDOFEquationNumbers();

  /**
   * Gather the DOF value on the root proccessor
   *
   * @param to_gather data to gather
   * @param root processor on which data are gathered
   * @param gathered Vector containing the gathered data, only valid on root processor
   */
  template<typename T>
  void gather(const Vector<T> & to_gather, UInt root,
	      Vector<T> * gathered = NULL) const;

  /**
   * Scatter a DOF Vector form root to all processors
   *
   * @param scattered data to scatter, only valid on root processor
   * @param root processor scattering data
   * @param to_scatter result of scattered data
   */
  template<typename T>
  void scatter(Vector<T> & scattered, UInt root,
	       const Vector<T> * to_scatter = NULL) const;


  template<typename T> void synchronize(Vector<T> & vector) const ;
  template<class Op, typename T> void reduceSynchronize(Vector<T> & vector) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the equation_number Vector
  AKANTU_GET_MACRO(LocalDOFEquationNumbers, local_dof_equation_numbers, const Vector<Int> &);
  AKANTU_GET_MACRO(GlobalDOFEquationNumbers, global_dof_equation_numbers, const Vector<Int> &);

  Vector<Int> * getLocalDOFEquationNumbersPointer(){return &local_dof_equation_numbers;};
  Vector<Int> * getGlobalDOFEquationNumbersPointer(){return &global_dof_equation_numbers;};

  typedef unordered_map<Int, UInt>::type GlobalEquationNumberMap;
  AKANTU_GET_MACRO(GlobalEquationNumberToLocal, global_dof_equation_number_to_local, const GlobalEquationNumberMap &)

  /// get the Vector of global ids of the dofs
  AKANTU_GET_MACRO(DOFGlobalIDs, dof_global_ids, const Vector<UInt> &);

  /// get the global id of a dof
  inline UInt getDOFGlobalID(UInt local_id) const { return dof_global_ids(local_id); }

  /// get the local id of a global dof
  inline UInt getDOFLocalID(UInt global_id) const { return global_dof_to_local.find(global_id)->second; }

  /// get the DOF type Vector
  AKANTU_GET_MACRO(DOFTypes, dof_types, const Vector<Int> &);

  AKANTU_GET_MACRO(NbDOFs, nb_dofs, UInt);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// equation number position where a dof is synchronized in the matrix (by default = global id)
  Vector<Int> global_dof_equation_numbers;
  Vector<Int> local_dof_equation_numbers;
  GlobalEquationNumberMap global_dof_equation_number_to_local;

  /// DOF global id
  Vector<UInt> dof_global_ids;

  /*
   * DOF type  -3 pure ghost, -2  master for the dof, -1 normal dof,  i in
   * [0-N] slave dof and master is proc i
   */
  Vector<Int> dof_types;

  /// number of dofs
  UInt nb_dofs;

  UInt nb_global_dofs;

  unordered_map<UInt, UInt>::type global_dof_to_local;

  UInt prank;
  UInt psize;
  StaticCommunicator * communicator;

  struct PerProcInformations {
    /// dofs to send to the proc
    Vector<UInt> slave_dofs;
    /// dofs to recvs from the proc
    Vector<UInt> master_dofs;

    /* ---------------------------------------------------------------------- */
    /* Data for gather/scatter                                                */
    /* ---------------------------------------------------------------------- */
    /// the dof that the node handle
    Vector<UInt> dofs;
    /// the dof that the proc need
    Vector<UInt> needed_dofs;
  };

  std::vector<PerProcInformations> proc_informations;

  /// nb dofs with type -1 or -2
  UInt nb_local_dofs;
  /// nb dof with type >= 0
  UInt nb_needed_dofs;

  bool gather_scatter_scheme_initialized;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "dof_synchronizer_inline_impl.cc"
#endif

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const DOFSynchronizer & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__

#endif /* __AKANTU_DOF_SYNCHRONIZER_HH__ */
