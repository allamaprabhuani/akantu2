/**
 * @file   material.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:43:41 2010
 *
 * @brief  Implementation of the common part of the material class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "material.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialBase::MaterialBase(SolidMechanicsModel & model, const MaterialID & id) :
  Memory(model.getMemoryID()), id(id), model(&model) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
