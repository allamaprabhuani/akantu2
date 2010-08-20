/**
 * @file   solid_mechanics_model_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 29 12:07:04 2010
 *
 * @brief  Implementation of the inline functions of the SolidMechanicsModel class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
// inline Vector<Real> & SolidMechanicsModel::getStress(ElementType type) {
//   AKANTU_DEBUG_IN();
//   AKANTU_DEBUG_ASSERT(stress[type] != NULL,
// 		      "The model " << id << " has no element of kind : "<< type);

//   AKANTU_DEBUG_OUT();
//   return *(stress[type]);
// }

// /* -------------------------------------------------------------------------- */
// inline Vector<Real> & SolidMechanicsModel::getStrain(ElementType type) {
//   AKANTU_DEBUG_IN();
//   AKANTU_DEBUG_ASSERT(strain[type] != NULL,
// 		      "The model " << id << " has no element of kind : "<< type);
//   AKANTU_DEBUG_OUT();
//   return *(strain[type]);
// }

/* -------------------------------------------------------------------------- */
inline Material & SolidMechanicsModel::getMaterial(UInt mat_index) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(mat_index < materials.size(),
		      "The model " << id << " has no material no "<< mat_index);
  AKANTU_DEBUG_OUT();
  return *materials[mat_index];
}
