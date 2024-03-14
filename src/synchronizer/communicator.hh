/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "aka_common.hh"
#include "aka_event_handler_manager.hh"
#include "communication_buffer.hh"
#include "communication_request.hh"
#include "communicator_event_handler.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_STATIC_COMMUNICATOR_HH_
#define AKANTU_STATIC_COMMUNICATOR_HH_

namespace akantu {

namespace debug {
  class CommunicationException : public Exception {
  public:
    CommunicationException()
        : Exception("An exception happen during a communication process.") {}
  };
} // namespace debug

/// @enum SynchronizerOperation reduce operation that the synchronizer can
/// perform
enum class SynchronizerOperation {
  _sum,
  _min,
  _max,
  _prod,
  _land,
  _band,
  _lor,
  _bor,
  _lxor,
  _bxor,
  _min_loc,
  _max_loc,
  _null
};

enum class CommunicationMode { _auto, _synchronous, _ready };

namespace {
  constexpr int _any_source = -1;
}
} // namespace akantu

namespace akantu {

struct CommunicatorInternalData {
  CommunicatorInternalData() = default;
  CommunicatorInternalData(const CommunicatorInternalData &) = default;
  CommunicatorInternalData(CommunicatorInternalData &&) = default;
  CommunicatorInternalData &
  operator=(const CommunicatorInternalData &) = default;
  CommunicatorInternalData & operator=(CommunicatorInternalData &&) = default;
  virtual ~CommunicatorInternalData() = default;
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

class Communicator : public EventHandlerManager<CommunicatorEventHandler> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  struct private_member {
    private_member() = default;
    private_member(const private_member &) = default;
    private_member(private_member &&) = default;
    private_member & operator=(const private_member &) = default;
    private_member & operator=(private_member &&) = default;
    virtual ~private_member() = default;
  };

  [[deprecated("argc, argv arguments are not nedded anymore")]] Communicator(
      int & argc, char **& argv, const private_member & /*m*/);
  Communicator(const private_member & /*unused*/ = private_member{});
  ~Communicator() override;
  Communicator(const Communicator &) = delete;
  Communicator(Communicator &&) = delete;
  Communicator & operator=(const Communicator &) = delete;
  Communicator & operator=(Communicator &&) = delete;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Point to Point                                                           */
  /* ------------------------------------------------------------------------ */
  template <typename T>
  void probe(Int sender, Int tag, CommunicationStatus & status) const;

  template <typename T>
  [[nodiscard]] bool asyncProbe(Int sender, Int tag,
                                CommunicationStatus & status) const;

  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void receive(Array<T> & values, Int sender, Int tag) const {
    return this->receiveImpl(
        values.data(), values.size() * values.getNbComponent(), sender, tag);
  }

  template <typename T>
  inline void receive(std::vector<T> & values, Int sender, Int tag) const {
    return this->receiveImpl(values.data(), values.size(), sender, tag);
  }

  template <typename Tensor>
  inline void
  receive(Tensor & values, Int sender, Int tag,
          std::enable_if_t<aka::is_tensor<Tensor>::value> * /*unused*/ =
              nullptr) const {
    return this->receiveImpl(values.data(), values.size(), sender, tag);
  }

  inline void receive(CommunicationBufferTemplated<true> & values, Int sender,
                      Int tag) const {
    return this->receiveImpl(values.data(), Int(values.size()), sender, tag);
  }

  inline void receive(CommunicationBufferTemplated<false> & values, Int sender,
                      Int tag) const {
    CommunicationStatus status;
    this->probe<char>(sender, tag, status);
    values.reserve(status.size());
    return this->receiveImpl(values.data(), int(values.size()), sender, tag);
  }

  template <typename T>
  inline void
  receive(T & values, Int sender, Int tag,
          std::enable_if_t<std::is_arithmetic<T>::value> * /*unused*/ =
              nullptr) const {
    return this->receiveImpl(&values, 1, sender, tag);
  }
  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void
  send(const Array<T> & values, Int receiver, Int tag,
       const CommunicationMode & mode = CommunicationMode::_auto) const {
    return this->sendImpl(values.data(),
                          values.size() * values.getNbComponent(), receiver,
                          tag, mode);
  }

  template <typename T>
  inline void
  send(const std::vector<T> & values, Int receiver, Int tag,
       const CommunicationMode & mode = CommunicationMode::_auto) const {
    return this->sendImpl(values.data(), values.size(), receiver, tag, mode);
  }

  template <typename Tensor,
            std::enable_if_t<aka::is_tensor_v<Tensor>> * = nullptr>
  inline void
  send(const Tensor & values, Int receiver, Int tag,
       const CommunicationMode & mode = CommunicationMode::_auto) const {
    return this->sendImpl(values.data(), values.size(), receiver, tag, mode);
  }

  template <bool is_static>
  inline void
  send(const CommunicationBufferTemplated<is_static> & values, Int receiver,
       Int tag,
       const CommunicationMode & mode = CommunicationMode::_auto) const {
    return this->sendImpl(values.data(), values.size(), receiver, tag, mode);
  }
  template <typename T>
  inline void send(const T & values, Int receiver, Int tag,
                   const CommunicationMode & mode = CommunicationMode::_auto,
                   std::enable_if_t<std::is_arithmetic<T>::value> * /*unused*/ =
                       nullptr) const {
    return this->sendImpl(&values, 1, receiver, tag, mode);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  [[nodiscard]] inline CommunicationRequest
  asyncSend(const Array<T> & values, Int receiver, Int tag,
            const CommunicationMode & mode = CommunicationMode::_auto) const {
    return this->asyncSendImpl(values.data(),
                               values.size() * values.getNbComponent(),
                               receiver, tag, mode);
  }
  template <typename T>
  [[nodiscard]] inline CommunicationRequest
  asyncSend(const std::vector<T> & values, Int receiver, Int tag,
            const CommunicationMode & mode = CommunicationMode::_auto) const {
    return this->asyncSendImpl(values.data(), values.size(), receiver, tag,
                               mode);
  }

  template <typename Tensor>
  [[nodiscard]] inline CommunicationRequest asyncSend(
      const Tensor & values, Int receiver, Int tag,
      const CommunicationMode & mode = CommunicationMode::_auto,
      std::enable_if_t<aka::is_tensor_v<Tensor>> * /*unused*/ = nullptr) const {
    return this->asyncSendImpl(values.data(), values.size(), receiver, tag,
                               mode);
  }
  template <bool is_static>
  [[nodiscard]] inline CommunicationRequest
  asyncSend(const CommunicationBufferTemplated<is_static> & values,
            Int receiver, Int tag,
            const CommunicationMode & mode = CommunicationMode::_auto) const {
    return this->asyncSendImpl(values.data(), values.size(), receiver, tag,
                               mode);
  }
  template <typename T>
  [[nodiscard]] inline CommunicationRequest
  asyncSend(const T & values, Int receiver, Int tag,
            const CommunicationMode & mode = CommunicationMode::_auto,
            std::enable_if_t<std::is_arithmetic<T>::value> * /*unused*/ =
                nullptr) const {
    return this->asyncSendImpl(&values, 1, receiver, tag, mode);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  [[nodiscard]] inline CommunicationRequest
  asyncReceive(Array<T> & values, Int sender, Int tag) const {
    return this->asyncReceiveImpl(
        values.data(), values.size() * values.getNbComponent(), sender, tag);
  }
  template <typename T>
  [[nodiscard]] inline CommunicationRequest
  asyncReceive(std::vector<T> & values, Int sender, Int tag) const {
    return this->asyncReceiveImpl(values.data(), values.size(), sender, tag);
  }

  template <typename Tensor,
            typename = std::enable_if_t<aka::is_tensor_v<Tensor>>>
  [[nodiscard]] inline CommunicationRequest
  asyncReceive(Tensor & values, Int sender, Int tag) const {
    return this->asyncReceiveImpl(values.data(), values.size(), sender, tag);
  }
  template <bool is_static>
  [[nodiscard]] inline CommunicationRequest
  asyncReceive(CommunicationBufferTemplated<is_static> & values, Int sender,
               Int tag) const {
    return this->asyncReceiveImpl(values.data(), values.size(), sender, tag);
  }

  /* ------------------------------------------------------------------------ */
  /* Collectives                                                              */
  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void
  allReduce(Array<T> & values,
            SynchronizerOperation op = SynchronizerOperation::_sum) const {
    this->allReduceImpl(values.data(), values.size() * values.getNbComponent(),
                        op);
  }

  template <typename T>
  inline void
  allReduce(std::vector<T> & values,
            SynchronizerOperation op = SynchronizerOperation::_sum) const {
    this->allReduceImpl(values.data(), values.size(), op);
  }

  template <typename Derived>
  inline void
  allReduce(Eigen::MatrixBase<Derived> & values,
            SynchronizerOperation op = SynchronizerOperation::_sum) const {
    this->allReduceImpl(values.derived().data(), values.size(), op);
  }

  template <typename T>
  inline void
  allReduce(T & values, SynchronizerOperation op = SynchronizerOperation::_sum,
            std::enable_if_t<std::is_arithmetic<T>::value> * /*unused*/ =
                nullptr) const {
    this->allReduceImpl(&values, 1, op);
  }

  template <typename T>
  inline void
  scan(Array<T> & values,
       SynchronizerOperation op = SynchronizerOperation::_sum) const {
    this->scanImpl(values.data(), values.size() * values.getNbComponent(), op);
  }

  template <typename Tensor>
  inline void scan(
      Tensor & values, SynchronizerOperation op,
      std::enable_if_t<aka::is_tensor_v<Tensor>> * /*unused*/ = nullptr) const {
    this->scanImpl(values.data(), values.data(), values.size(), op);
  }

  template <typename T>
  inline void scan(T & values,
                   SynchronizerOperation op = SynchronizerOperation::_sum,
                   std::enable_if_t<std::is_arithmetic<T>::value> * /*unused*/ =
                       nullptr) const {
    this->scanImpl(&values, &values, 1, op);
  }

  template <typename T>
  inline void
  exclusiveScan(Array<T> & values,
                SynchronizerOperation op = SynchronizerOperation::_sum) const {
    this->exclusiveScanImpl(values.data(), values.data(),
                            values.size() * values.getNbComponent(), op);
  }

  template <typename Tensor,
            std::enable_if_t<aka::is_tensor_v<Tensor>> * = nullptr>
  inline void
  exclusiveScan(Tensor & values,
                SynchronizerOperation op = SynchronizerOperation::_sum) const {
    this->exclusiveScanImpl(values.data(), values.size(), op);
  }

  template <typename T, size_t size>
  inline void
  exclusiveScan(std::array<T, size> & values,
                SynchronizerOperation op = SynchronizerOperation::_sum) const {
    this->exclusiveScanImpl(values.data(), values.data(), size, op);
  }

  template <typename T>
  inline void
  exclusiveScan(std::vector<T> & values, std::vector<T> & results,
                SynchronizerOperation op = SynchronizerOperation::_sum) const {
    results.resize(values.size());
    this->exclusiveScanImpl(values.data(), results.data(), values.size(), op);
  }

  template <typename T,
            std::enable_if_t<std::is_arithmetic<T>::value> * = nullptr>
  inline void
  exclusiveScan(T & values,
                SynchronizerOperation op = SynchronizerOperation::_sum) const {
    this->exclusiveScanImpl(&values, &values, 1, op);
  }

  template <typename T,
            std::enable_if_t<std::is_arithmetic<T>::value> * = nullptr>
  inline void
  exclusiveScan(T & values, T & result,
                SynchronizerOperation op = SynchronizerOperation::_sum) const {
    this->exclusiveScanImpl(&values, &result, 1, op);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T> inline void allGather(Array<T> & values) const {
    AKANTU_DEBUG_ASSERT(getNbProc() == values.size(),
                        "The array size is not correct");
    this->allGatherImpl(values.data(), values.getNbComponent());
  }

  template <typename Tensor,
            typename = std::enable_if_t<aka::is_tensor_v<Tensor>>>
  inline void allGather(Tensor & values) const {
    AKANTU_DEBUG_ASSERT(values.size() / getNbProc() > 0,
                        "The vector size is not correct");
    this->allGatherImpl(values.data(), values.size() / getNbProc());
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void allGatherV(Array<T> & values, const Array<Int> & sizes) const {
    this->allGatherVImpl(values.data(), sizes.data());
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void reduce(Array<T> & values, SynchronizerOperation op,
                     int root = 0) const {
    this->reduceImpl(values.data(), values.size() * values.getNbComponent(), op,
                     root);
  }

  /* ------------------------------------------------------------------------ */
  template <typename Tensor>
  inline void gather(
      Tensor & values, int root = 0,
      std::enable_if_t<aka::is_tensor_v<Tensor>> * /*unused*/ = nullptr) const {
    this->gatherImpl(values.data(), values.getNbComponent(), root);
  }
  template <typename T>
  inline void
  gather(T values, int root = 0,
         std::enable_if_t<std::is_arithmetic<T>::value> * /*unused*/ =
             nullptr) const {
    this->gatherImpl(&values, 1, root);
  }
  /* ------------------------------------------------------------------------ */
  template <typename Tensor, typename T>
  inline void gather(
      Tensor & values, Array<T> & gathered,
      std::enable_if_t<aka::is_tensor_v<Tensor>> * /*unused*/ = nullptr) const {
    AKANTU_DEBUG_ASSERT(values.size() == gathered.getNbComponent(),
                        "The array size is not correct");
    gathered.resize(getNbProc());
    this->gatherImpl(values.data(), values.size(), gathered.data(),
                     gathered.getNbComponent());
  }

  template <typename T>
  inline void
  gather(T values, Array<T> & gathered,
         std::enable_if_t<std::is_arithmetic<T>::value> * /*unused*/ =
             nullptr) const {
    this->gatherImpl(&values, 1, gathered.data(), 1);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void gatherV(Array<T> & values, const Array<Int> & sizes,
                      int root = 0) const {
    this->gatherVImpl(values.data(), sizes.data(), root);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void broadcast(Array<T> & values, int root = 0) const {
    this->broadcastImpl(values.data(),
                        int(values.size() * values.getNbComponent()), root);
  }

  template <typename T>
  inline void broadcast(std::vector<T> & values, int root = 0) const {
    this->broadcastImpl(values.data(), int(values.size()), root);
  }

  inline void broadcast(CommunicationBufferTemplated<true> & buffer,
                        int root = 0) const {
    this->broadcastImpl(buffer.data(), int(buffer.size()), root);
  }

  inline void broadcast(CommunicationBufferTemplated<false> & buffer,
                        int root = 0) const {
    auto buffer_size = buffer.size();
    this->broadcastImpl(&buffer_size, 1, root);
    if (whoAmI() != root) {
      buffer.reserve(buffer_size);
    }

    if (buffer_size == 0) {
      return;
    }
    this->broadcastImpl(buffer.data(), int(buffer.size()), root);
  }

  template <typename T> inline void broadcast(T & values, int root = 0) const {
    this->broadcastImpl(&values, 1, root);
  }

  /* ------------------------------------------------------------------------ */
  void barrier() const;
  [[nodiscard]] CommunicationRequest asyncBarrier() const;

  /* ------------------------------------------------------------------------ */
  /* Request handling                                                         */
  /* ------------------------------------------------------------------------ */
  static bool test(CommunicationRequest & request);
  static bool testAll(std::vector<CommunicationRequest> & request);
  static void wait(CommunicationRequest & request);
  static void waitAll(std::vector<CommunicationRequest> & requests);
  static Int waitAny(std::vector<CommunicationRequest> & requests);
  static inline void freeCommunicationRequest(CommunicationRequest & request);
  static inline void
  freeCommunicationRequest(std::vector<CommunicationRequest> & requests);

  template <typename T, typename MsgProcessor>
  inline void
  receiveAnyNumber(std::vector<CommunicationRequest> & send_requests,
                   MsgProcessor && processor, Int tag) const;

protected:
  template <typename T>
  void
  sendImpl(const T * buffer, Int size, Int receiver, Int tag,
           const CommunicationMode & mode = CommunicationMode::_auto) const;

  template <typename T>
  void receiveImpl(T * buffer, Int size, Int sender, Int tag) const;

  template <typename T>
  CommunicationRequest asyncSendImpl(
      const T * buffer, Int size, Int receiver, Int tag,
      const CommunicationMode & mode = CommunicationMode::_auto) const;

  template <typename T>
  CommunicationRequest asyncReceiveImpl(T * buffer, Int size, Int sender,
                                        Int tag) const;

  template <typename T>
  void allReduceImpl(T * values, int nb_values, SynchronizerOperation op) const;

  template <typename T>
  void scanImpl(T * values, T * results, int nb_values,
                SynchronizerOperation op) const;

  template <typename T>
  void exclusiveScanImpl(T * values, T * results, int nb_values,
                         SynchronizerOperation op) const;

  template <typename T> void allGatherImpl(T * values, int nb_values) const;
  template <typename T>
  void allGatherVImpl(T * values, const int * nb_values) const;

  template <typename T>
  void reduceImpl(T * values, int nb_values, SynchronizerOperation op,
                  int root = 0) const;
  template <typename T>
  void gatherImpl(T * values, int nb_values, int root = 0) const;

  template <typename T>
  void gatherImpl(T * values, int nb_values, T * gathered,
                  int nb_gathered = 0) const;

  template <typename T>
  void gatherVImpl(T * values, int * nb_values, int root = 0) const;

  template <typename T>
  void broadcastImpl(T * values, int nb_values, int root = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  [[nodiscard]] Int getNbProc() const;
  [[nodiscard]] Int whoAmI() const;

  [[deprecated("Use getWorldCommunicator or getSelfCommunicator "
               "instead")]] static Communicator &
  getStaticCommunicator();
  [[deprecated("Use getWorldCommunicator or getSelfCommunicator "
               "instead")]] static Communicator &
  getStaticCommunicator(int & argc, char **& argv);

  static Communicator & getSelfCommunicator();
  static Communicator & getWorldCommunicator();
  static Communicator & getNullCommunicator();

  [[nodiscard]] int getMaxTag() const;
  [[nodiscard]] int getMinTag() const;

  AKANTU_GET_MACRO(CommunicatorData, (*communicator_data), decltype(auto));

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  static std::unique_ptr<Communicator>
      world_communicator; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)
  static std::unique_ptr<Communicator>
      self_communicator; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)
  static std::unique_ptr<Communicator>
      null_communicator; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)

protected:
  std::unique_ptr<CommunicatorInternalData> communicator_data;
};

inline std::ostream & operator<<(std::ostream & stream,
                                 const CommunicationRequest & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "communicator_inline_impl.hh"

#endif /* AKANTU_STATIC_COMMUNICATOR_HH_ */
