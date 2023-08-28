#include "aka_common.hh"

#ifndef AKANTU_RELEASE_HH_
#define AKANTU_RELEASE_HH_

namespace akantu {

struct Release {
  Release() = default;
  Release(Int release) : release(release) {}

  Release & operator--() {
    --release;
    return *this;
  }

  Release & operator++() {
    ++release;
    if (release == -1) { // -1 is the special state always different
      ++release;
    }
    return *this;
  }

  Release operator++(int) {
    auto tmp = *this;
    ++(*this);
    return tmp;
  }

  bool operator!=(const Release & other) const {
    return (release == -1) or (other.release == -1) or
           (other.release != release);
  }

  bool operator==(const Release & other) const { return not(*this != other); }

  operator Int() const { return release; }

private:
  Int release{-1};
};

} // namespace akantu

#endif // AKANTU_RELEASE_HH_
