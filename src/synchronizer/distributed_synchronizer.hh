/**
 * @file   distributed_synchronizer.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Jun 16 16:36:52 2011
 *
 * @brief  wrapper to the static communicator
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

#ifndef __AKANTU_DISTRIBUTED_SYNCHRONIZER_HH__
#define __AKANTU_DISTRIBUTED_SYNCHRONIZER_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_vector.hh"
#include "static_communicator.hh"
#include "synchronizer.hh"
#include "mesh.hh"
#include "mesh_partition.hh"
#include "communication_buffer.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class DistributedSynchronizer : public Synchronizer, public MeshEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DistributedSynchronizer(Mesh & mesh,
                          SynchronizerID id = "distributed_synchronizer",
                          MemoryID memory_id = 0);

public:
  virtual ~DistributedSynchronizer();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// get a  mesh and a partition and  create the local mesh  and the associated
  /// DistributedSynchronizer
  static DistributedSynchronizer *
  createDistributedSynchronizerMesh(Mesh & mesh,
                                    const MeshPartition * partition,
                                    UInt root = 0,
                                    SynchronizerID id = "distributed_synchronizer",
                                    MemoryID memory_id = 0);

  /* ------------------------------------------------------------------------ */
  /* Inherited from Synchronizer                                              */
  /* ------------------------------------------------------------------------ */

  /// asynchronous synchronization of ghosts
  void asynchronousSynchronize(DataAccessor & data_accessor,SynchronizationTag tag);

  /// wait end of asynchronous synchronization of ghosts
  void waitEndSynchronize(DataAccessor & data_accessor,SynchronizationTag tag);

  /// build processor to element corrispondance
  void buildPrankToElement(ByElementTypeUInt & prank_to_element);

  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// mesh event handler onRemovedElement
  virtual void onElementsRemoved(const Array<Element> & element_list,
                                 const ByElementTypeUInt & new_numbering,
                                 const RemovedElementsEvent & event);

protected:
  /// fill the nodes type vector
  void fillNodesType(Mesh & mesh);

void fillNodesType(const MeshData & mesh_data,
                   DynamicCommunicationBuffer * buffers,
                   const std::string & tag_name,
                   const ElementType & el_type,
                   const UInt * partition_num);

template<typename T>
void fillTagBufferTemplated(const MeshData & mesh_data,
                            DynamicCommunicationBuffer * buffers,
                            const std::string & tag_name,
                            const ElementType & el_type,
                            const UInt * partition_num,
                            const UInt * ghost_partition,
                            const UInt * ghost_partition_offset);

void fillTagBuffer(const MeshData & mesh_data,
                   DynamicCommunicationBuffer * buffers,
                   const std::string & tag_name,
                   const ElementType & el_type,
                   const UInt * partition_num,
                   const UInt * ghost_partition,
                   const UInt * ghost_partition_offset);

template<typename T, typename BufferType>
void populateMeshDataTemplated(MeshData & mesh_data,
                               BufferType & buffer,
                               const std::string & tag_name,
                               const ElementType & el_type,
                               UInt nb_component,
                               UInt nb_local_element,
                               UInt nb_ghost_element);

template <typename BufferType>
void populateMeshData(MeshData & mesh_data,
                      BufferType & buffer,
                      const std::string & tag_name,
                      const ElementType & el_type,
                      const MeshDataTypeCode & type_code,
                      UInt nb_component,
                      UInt nb_local_element,
                      UInt nb_ghost_element);

  /// fill the communications array of a distributedSynchronizer based on a partition array
  void fillCommunicationScheme(UInt * partition,
                               UInt nb_local_element,
                               UInt nb_ghost_element,
                               ElementType type);

  /// compute buffer size for a given tag and data accessor
  void computeBufferSize(DataAccessor & data_accessor, SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  enum CommTags {
    TAG_SIZES        = 0,
    TAG_CONNECTIVITY = 1,
    TAG_DATA         = 2,
    TAG_PARTITIONS   = 3,
    TAG_NB_NODES     = 4,
    TAG_NODES        = 5,
    TAG_COORDINATES  = 6,
    TAG_NODES_TYPE   = 7,
    TAG_MESH_DATA   = 8
  };

protected:
  /// reference to the underlying mesh
  Mesh & mesh;

  /// the static memory instance
  StaticCommunicator * static_communicator;

  class Communication {
  public:
    void resize(UInt size) {
      send_buffer.resize(size);
      recv_buffer.resize(size);
      size_to_send   .resize(size);
      size_to_receive.resize(size);
    }

  public:
    /// size of data to send to each processor
    std::vector<UInt> size_to_send;
    /// size of data to recv to each processor
    std::vector<UInt> size_to_receive;
    std::vector< CommunicationBuffer > send_buffer;
    std::vector< CommunicationBuffer > recv_buffer;

    std::vector<CommunicationRequest *> send_requests;
    std::vector<CommunicationRequest *> recv_requests;
  };

  std::map<SynchronizationTag, Communication> communications;

  /// list of element to send to proc p
  Array<Element> * send_element;
  /// list of element to receive from proc p
  Array<Element> * recv_element;

  UInt nb_proc;
  UInt rank;

  friend class FacetSynchronizer;

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
template<typename T>
void DistributedSynchronizer::fillTagBufferTemplated(const MeshData & mesh_data,
                                                     DynamicCommunicationBuffer * buffers,
                                                     const std::string & tag_name,
                                                     const ElementType & el_type,
                                                     const UInt * partition_num,
                                                     const UInt * ghost_partition,
                                                     const UInt * ghost_partition_offset) {
  const Array<T> & data = mesh_data.getElementalDataArray<T>(tag_name, el_type);
  // Not possible to use the iterator because it potentially triggers the creation of complex
  // type templates (such as akantu::Vector< std::vector<Element> > which don't implement the right interface
  // (e.g. operator<< in that case).
  //typename Array<T>::template const_iterator< Vector<T> > data_it  = data.begin(data.getNbComponent());
  //typename Array<T>::template const_iterator< Vector<T> > data_end = data.end(data.getNbComponent());

  const T * data_it = data.storage();
  const T * data_end = data.storage() + data.getSize()*data.getNbComponent();
  const UInt * part = partition_num;

  /// copying the data, element by element
  for (; data_it != data_end; ++part) {
    for(UInt j(0); j < data.getNbComponent(); ++j, ++data_it) {
      buffers[*part] << *data_it;
    }
  }

  data_it  = data.storage();
  const UInt * offset = ghost_partition_offset;
  /// copying the data for the ghost element
  for (; data_it != data_end; data_it+=data.getNbComponent(), ++offset) {
    for (UInt p = *offset; p < *(offset + 1); ++p) {
      UInt proc = ghost_partition[p];
      for(UInt j(0); j < data.getNbComponent(); ++j) {
        buffers[proc] << data_it[j];
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename BufferType>
void DistributedSynchronizer::populateMeshData(MeshData & mesh_data,
                                               BufferType & buffer,
                                               const std::string & tag_name,
                                               const ElementType & el_type,
                                               const MeshDataTypeCode & type_code,
                                               UInt nb_component,
                                               UInt nb_local_element,
                                               UInt nb_ghost_element) {
  #define AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA(r, extra_param, elem)	\
    case BOOST_PP_TUPLE_ELEM(2, 0, elem) : { \
      populateMeshDataTemplated<BOOST_PP_TUPLE_ELEM(2, 1, elem)>(mesh_data, buffer, tag_name, el_type, nb_component, nb_local_element, nb_ghost_element); \
      break; \
    } \

  switch(type_code) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA, , AKANTU_MESH_DATA_TYPES)
  default : AKANTU_DEBUG_ERROR("Could not determine the type of tag" << tag_name << "!"); break;
  }
  #undef AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA
}

/* -------------------------------------------------------------------------- */
template<typename T, typename BufferType>
void DistributedSynchronizer::populateMeshDataTemplated(MeshData & mesh_data,
                                                        BufferType & buffer,
                                                        const std::string & tag_name,
                                                        const ElementType & el_type,
                                                        UInt nb_component,
                                                        UInt nb_local_element,
                                                        UInt nb_ghost_element) {

  if(nb_local_element != 0) {
    mesh_data.registerElementalData<T>(tag_name);
    Array<T> & data = mesh_data.getElementalDataArrayAlloc<T>(tag_name, el_type, _not_ghost, nb_component);
    data.resize(nb_local_element);
    /// unpacking the data, element by element
    for (UInt i(0); i < nb_local_element; ++i) {
      for(UInt j(0); j < nb_component; ++j) {
        buffer >> data(i,j);
      }
    }
  }

  if(nb_ghost_element != 0) {
    mesh_data.registerElementalData<T>(tag_name);
    Array<T> & data_ghost = mesh_data.getElementalDataArrayAlloc<T>(tag_name, el_type, _ghost, nb_component);
    data_ghost.resize(nb_ghost_element);

    /// unpacking the ghost data, element by element
    for (UInt j(0); j < nb_ghost_element; ++j) {
      for(UInt k(0); k < nb_component; ++k) {
        buffer >> data_ghost(j, k);
      }
    }
  }
}

//#include "distributedSynchronizer_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_DISTRIBUTED_SYNCHRONIZER_HH__ */
