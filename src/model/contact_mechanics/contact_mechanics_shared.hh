#ifndef __AKANTU_CONTACT_MECHANICS_SHARED_HH__
#define __AKANTU_CONTACT_MECHANICS_SHARED_HH__

#include "aka_common.hh"
#include "aka_grid_dynamic.hh"
#include "fe_engine.hh"
#include "model.hh"

/* -------------------------------------------------------------------------- */
/* Shared abstractions across different contact mechanics models              */
/* -------------------------------------------------------------------------- */

namespace akantu {

/// Base class for contact detectors
class AbstractContactDetector {

  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructors                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  AbstractContactDetector(Mesh & mesh, Array<Real> initial_positions);

public:
  virtual ~AbstractContactDetector() = default;

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// Return a new vector with the positions of a node.
  inline Vector<Real> getNodePosition(UInt node) const;

protected:
  /// Fill the matrix with the coordinates of an element.
  inline void coordinatesOfElement(const Element & el, Matrix<Real> & coords) const;

  /// Compute the minimum and maximum element sizes.
  /// Make sure to call fillNodesToElements on the mesh before.
  template <typename vector_type>
  inline auto computeElementSizes(vector_type && nodes) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the mesh
  AKANTU_GET_MACRO(Mesh, mesh, Mesh &)

  AKANTU_GET_MACRO_NOT_CONST(Positions, positions, Array<Real> &);
  AKANTU_SET_MACRO(Positions, positions, Array<Real>);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Mesh
  Mesh & mesh;

  /// dimension of the model
  const UInt spatial_dimension;

  /// contains the updated positions of the nodes
  Array<Real> positions;
};

/// Base class for contact mechanics models
// TODO: evaluate whether this is _really_ useful or not...
class AbstractContactMechanicsModel : public Model {
public:
  AbstractContactMechanicsModel(Mesh & mesh, const ModelType & type,
        std::shared_ptr<DOFManager> dof_manager, UInt dim, const ID id);

  /// Perform the contact search.
  /// This is called by the coupler before the relevant terms are assembled, to
  /// update the internal state of the contact model.
  virtual void search() = 0;
};

} // namespace akantu

#include "contact_mechanics_shared_inline_impl.cc"

#endif // __AKANTU_CONTACT_MECHANICS_SHARED_HH__
