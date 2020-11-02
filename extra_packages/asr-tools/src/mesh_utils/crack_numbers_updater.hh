/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_CRACK_NUMBERS_UPDATER_HH__
#define __AKANTU_CRACK_NUMBERS_UPDATER_HH__

/* -------------------------------------------------------------------------- */
#include "data_accessor.hh"
#include "solid_mechanics_model.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {
class ElementSynchronizer;
} // namespace akantu

namespace akantu {

class CrackNumbersUpdater : public DataAccessor<Element> {
public:
  CrackNumbersUpdater(SolidMechanicsModel & model,
                      const ElementSynchronizer & synchronizer)
      : model(model), synchronizer(synchronizer) {}

  void communicateCrackNumbers();

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
private:
  /// Reference to the model
  SolidMechanicsModel & model;

  /// distributed synchronizer to communicate nodes positions
  const ElementSynchronizer & synchronizer;
};

} // namespace akantu

#include "crack_numbers_updater_inline_impl.cc"

#endif /* __AKANTU_CRACK_NUMBERS_UPDATER_HH__ */
