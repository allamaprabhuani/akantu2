/**
 * @file   material_damage_iterative.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Specialization of the class material damage to damage only one gauss
 * point at a time and propagate damage in a linear way. Max principal stress
 * criterion is used as a failure criterion.
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_damage_iterative.hh"
#include "communicator.hh"
#include "data_accessor.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

INSTANTIATE_MATERIAL(damage_iterative, MaterialDamageIterative);

} // namespace akantu
