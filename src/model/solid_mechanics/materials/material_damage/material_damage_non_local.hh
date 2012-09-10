/**
 * @file   material_damage_non_local.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 16 18:28:00 2012
 *
 * @brief  interface for non local damage material
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

#ifndef __AKANTU_MATERIAL_DAMAGE_NON_LOCAL_HH__
#define __AKANTU_MATERIAL_DAMAGE_NON_LOCAL_HH__

__BEGIN_AKANTU__

template<UInt spatial_dimension,
	 template <UInt> class MaterialDamageLocal,
	 template <UInt> class WeightFunction = BaseWeightFunction>
class MaterialDamageNonLocal : public MaterialDamageLocal<spatial_dimension>,
			       public MaterialNonLocal<spatial_dimension, WeightFunction> {
public:
  typedef MaterialNonLocal<spatial_dimension, WeightFunction> MaterialNonLocalParent;
  typedef MaterialDamageLocal<spatial_dimension> MaterialDamageParent;

  MaterialDamageNonLocal(SolidMechanicsModel & model, const ID & id)  :
    Material(model, id),
    MaterialDamageParent(model, id), MaterialNonLocalParent(model, id) { };

  /* ------------------------------------------------------------------------ */
  virtual void initMaterial() {
    MaterialDamageParent::initMaterial();
    MaterialNonLocalParent::initMaterial();
  }

protected:
  /* -------------------------------------------------------------------------- */
  virtual void computeNonLocalStress(ElementType type, GhostType ghost_type = _not_ghost) = 0;

  /* ------------------------------------------------------------------------ */
  void computeNonLocalStress(GhostType ghost_type) {
    AKANTU_DEBUG_IN();
    Mesh::type_iterator it = this->model->getFEM().getMesh().firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last_type = this->model->getFEM().getMesh().lastType(spatial_dimension, ghost_type);
    for(; it != last_type; ++it) {
      computeNonLocalStress(*it, ghost_type);
    }
    this->updateDissipatedEnergy(ghost_type);
    AKANTU_DEBUG_OUT();
  }

public:
  /* ------------------------------------------------------------------------ */
  virtual bool parseParam(const std::string & key, const std::string & value,
			const ID & id) {
    return MaterialNonLocalParent::parseParam(key, value, id) ||
      MaterialDamageParent::parseParam(key, value, id);
  }

public:
  /* ------------------------------------------------------------------------ */
  virtual inline UInt getNbDataToPack(const Element & element,
				      SynchronizationTag tag) const {
    return MaterialNonLocalParent::getNbDataToPack(element, tag) +
      MaterialDamageParent::getNbDataToPack(element, tag);
  }

  virtual inline UInt getNbDataToUnpack(const Element & element,
					SynchronizationTag tag) const {
    return MaterialNonLocalParent::getNbDataToUnpack(element, tag) +
      MaterialDamageParent::getNbDataToUnpack(element, tag);
  }


  virtual inline void packData(CommunicationBuffer & buffer,
			       const Element & element,
			       SynchronizationTag tag) const {
    MaterialNonLocalParent::packData(buffer, element, tag);
    MaterialDamageParent::packData(buffer, element, tag);
  }

  virtual inline void unpackData(CommunicationBuffer & buffer,
				 const Element & element,
				 SynchronizationTag tag) {
    MaterialNonLocalParent::unpackData(buffer, element, tag);
    MaterialDamageParent::unpackData(buffer, element, tag);
  }
};

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_DAMAGE_NON_LOCAL_HH__ */
