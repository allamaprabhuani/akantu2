#include "diffusion_law.hh"
#include "heat_transfer_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
DiffusionLaw::DiffusionLaw(DiffusionModel & model, const ID & id,
                           const ID & fe_engine_id)
    : Parent(model, id, model.getSpatialDimension(), _ek_regular, fe_engine_id),
      grad_u(registerInternal("∇u", model.getSpatialDimension(), fe_engine_id)),
      d_grad_u(
          registerInternal("D∇u", model.getSpatialDimension(), fe_engine_id)),
      diffusivity(registerInternal("diffusivity",
                                   model.getSpatialDimension() *
                                       model.getSpatialDimension(),
                                   fe_engine_id)) {}

/* -------------------------------------------------------------------------- */
void DiffusionLaw::computeGradU(ElementType type, GhostType ghost_type) {
  if (grad_u.getRelease(type, ghost_type) != -1 and
      grad_u.getRelease(type, ghost_type) ==
          getHandler().getDiffusionRelease()) {
    return;
  }

  auto & elem_filter = getElementFilter(type, ghost_type);
  this->getFEEngine().gradientOnIntegrationPoints(
      getHandler().getDiffusion(), grad_u(type, ghost_type), 1, type,
      ghost_type, elem_filter);

  grad_u.getRelease(type, ghost_type) = getHandler().getDiffusionRelease();
}
/* -------------------------------------------------------------------------- */
void DiffusionLaw::computeDiffusivityGradU(GhostType ghost_type) {
  auto dim = getHandler().getSpatialDimension();
  for (auto && type : getElementFilter().elementTypes(dim, ghost_type)) {
    this->computeGradU(type, ghost_type);
    this->computeDiffusivityGradUOnQuadPoints(type, ghost_type);
  }
}

/* -------------------------------------------------------------------------- */
void DiffusionLaw::assembleInternalFlow(GhostType ghost_type) {
  auto && model = this->getHandler();
  auto dim = model.getSpatialDimension();
  auto & fem = getFEEngine();

  for (auto && type : getElementFilter().elementTypes(dim, ghost_type)) {
    auto && elem_filter = getElementFilter(type, ghost_type);
    auto nb_element = elem_filter.size();
    if (nb_element == 0) {
      return;
    }

    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    auto & d_gradu_vect = d_grad_u(type, ghost_type);
    auto nb_quad_points = d_grad_u.size();

    Array<Real> bt_d_gu(nb_quad_points, nb_nodes_per_element);
    fem.computeBtD(d_gradu_vect, bt_d_gu, type, ghost_type);

    Array<Real> int_bt_d_gu(nb_element, nb_nodes_per_element);

    fem.integrate(bt_d_gu, int_bt_d_gu, nb_nodes_per_element, type, ghost_type,
                  elem_filter);

    model.getDOFManager().assembleElementalArrayLocalArray(
        int_bt_d_gu, model.getInternalFlow(), type, ghost_type, -1,
        elem_filter);
  }
}

/* -------------------------------------------------------------------------- */
void DiffusionLaw::assembleDiffusivityMatrix() {
  auto & fem = this->getFEEngine();
  auto & model = this->getHandler();

  for (auto && type : getElementFilter().elementTypes(spatial_dimension)) {
    auto && elem_filter = getElementFilter(type, _not_ghost);
    auto nb_element = elem_filter.size();
    if (nb_element == 0) {
      return;
    }
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type);

    auto bt_d_b = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points,
        nb_nodes_per_element * nb_nodes_per_element, "B^t*D*B");

    fem.computeBtDB(diffusivity(type), *bt_d_b, 2, type, _not_ghost,
                    elem_filter);

    /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
    auto K_e = std::make_unique<Array<Real>>(
        nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_e");

    fem.integrate(*bt_d_b, *K_e, nb_nodes_per_element * nb_nodes_per_element,
                  type, _not_ghost, elem_filter);

    model.getDOFManager().assembleElementalMatricesToMatrix(
        "K", model.getDOFName(), *K_e, type, _not_ghost, _symmetric,
        elem_filter);
  }
}

/* -------------------------------------------------------------------------- */
Int DiffusionLaw::getNbData(const Array<Element> & elements,
                            const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  Int size = 0;

  switch (tag) {
  case SynchronizationTag::_diffusion_gradient: {
    // temperature gradient
    size += getHandler().getNbIntegrationPoints(elements) * spatial_dimension *
            Int(sizeof(Real));
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
void DiffusionLaw::packData(CommunicationBuffer & buffer,
                            const Array<Element> & elements,
                            const SynchronizationTag & tag) const {
  switch (tag) {
  case SynchronizationTag::_diffusion_gradient: {
    packElementalDataHelper(this->grad_u, buffer, elements,
                            this->getFEEngine());
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
void DiffusionLaw::unpackData(CommunicationBuffer & buffer,
                              const Array<Element> & elements,
                              const SynchronizationTag & tag) {
  switch (tag) {
  case SynchronizationTag::_diffusion_gradient: {
    unpackElementalDataHelper(grad_u, buffer, elements, getFEEngine());
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
DiffusionFactory & DiffusionLaw::getFactory() {
  return DiffusionFactory::getInstance();
}

} // namespace akantu
