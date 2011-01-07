/**
 * @file   static_communicator.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 19 15:39:47 2010
 *
 * @brief  implementation of the common part of the static communicator
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
#include "static_communicator.hh"
#include "static_communicator_dummy.hh"

#ifdef AKANTU_USE_MPI
#  include "static_communicator_mpi.hh"
#endif

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
bool StaticCommunicator::is_instantiated = false;
StaticCommunicator * StaticCommunicator::static_communicator = NULL;

UInt CommunicationRequest::counter = 0;

/* -------------------------------------------------------------------------- */
CommunicationRequest::CommunicationRequest(UInt source, UInt dest) :
  source(source), destination(dest) {
  this->id = counter++;
}

/* -------------------------------------------------------------------------- */
CommunicationRequest::~CommunicationRequest() {

}

/* -------------------------------------------------------------------------- */
void CommunicationRequest::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "CommunicationRequest [" << std::endl;
  stream << space << " + id          : " << id << std::endl;
  stream << space << " + source      : " << source << std::endl;
  stream << space << " + destination : " << destination << std::endl;
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
StaticCommunicator * StaticCommunicator::getStaticCommunicator(CommunicatorType type) {
  AKANTU_DEBUG_IN();

  if(type == _communicator_mpi) {
#ifdef AKANTU_USE_MPI
    if (!static_communicator)
      AKANTU_DEBUG_ERROR("You must call getStaticCommunicator(argc, argv) to create a MPI communicator");
#else
    AKANTU_DEBUG_ERROR("You must recompile with the AKANTU_USE_MPI option ON to use a MPI communicator.");
#endif
  }

  if (!static_communicator)
    static_communicator = new StaticCommunicatorDummy();

  is_instantiated = true;

  AKANTU_DEBUG_OUT();
  return static_communicator;

}

/* -------------------------------------------------------------------------- */
StaticCommunicator * StaticCommunicator::getStaticCommunicator(__attribute__ ((unused)) int * argc,
							       __attribute__ ((unused)) char *** argv,
  							       CommunicatorType type) {
  if(type == _communicator_mpi) {
#ifdef AKANTU_USE_MPI
    if (!static_communicator)
      static_communicator = dynamic_cast<StaticCommunicator *>(new StaticCommunicatorMPI(argc, argv));
#else
    AKANTU_DEBUG_ERROR("You must recompile with the AKANTU_USE_MPI option ON to use a MPI communicator.");
#endif
  }

  return getStaticCommunicator(type);
}

__END_AKANTU__
