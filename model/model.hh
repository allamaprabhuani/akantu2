/**
 * @file   model.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 22 11:02:42 2010
 *
 * @brief  Interface of a model
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

#ifndef __AKANTU_MODEL_HH__
#define __AKANTU_MODEL_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_memory.hh"
#include "mesh.hh"
#include "fem.hh"
#include "mesh_utils.hh"
#include "synchronizer_registry.hh"
#include "distributed_synchronizer.hh"
#include "static_communicator.hh"
#include "mesh_partition.hh"
#include "dof_synchronizer.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class Model : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  inline Model(const ModelID & id = "model",
	       const MemoryID & memory_id = 0);

  inline virtual ~Model();

  typedef std::map<std::string, FEM *> FEMMap;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void initModel() = 0;

  /// create the synchronizer registry object
  void createSynchronizerRegistry(DataAccessor * data_accessor);

  /// create a parallel synchronizer and distribute the mesh
  Synchronizer & createParallelSynch(MeshPartition * partition,
				     DataAccessor * data_accessor);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /// return the object hadling equation numbers
  AKANTU_GET_MACRO(DOFSynchronizer, *dof_synchronizer, const DOFSynchronizer &)

  /// synchronize the boundary in case of parallel run
  virtual void synchronizeBoundaries() {};

  /// return the fem object associated with a provided name
  inline FEM & getFEM(std::string name = "") const;

  /// return the fem boundary object associated with a provided name
  virtual FEM & getFEMBoundary(std::string name = "");

  /// register a fem object associated with name
  template <typename FEMClass> inline void registerFEMObject(const std::string & name,
							     Mesh & mesh,
							     UInt spatial_dimension);
  /// unregister a fem object associated with name
  inline void unRegisterFEMObject(const std::string & name);

  /// return the synchronizer registry 
  SynchronizerRegistry & getSynchronizerRegistry(){return  *synch_registry;};

protected:
  /// return the fem object associated with a provided name
  template <typename FEMClass>
  inline FEMClass & getFEMClass(std::string name = "") const;

  /// return the fem boundary object associated with a provided name
  template <typename FEMClass>
  inline FEMClass & getFEMClassBoundary(std::string name = "");

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id
  ModelID id;

  /// the main fem object present in all  models
  FEMMap fems;

  /// the fem object present in all  models for boundaries
  FEMMap fems_boundary;

  /// default fem object
  std::string default_fem;

  /// synchronizer registry 
  SynchronizerRegistry * synch_registry;

  /// handle the equation number things
  DOFSynchronizer * dof_synchronizer;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "model_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const Model & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_MODEL_HH__ */
