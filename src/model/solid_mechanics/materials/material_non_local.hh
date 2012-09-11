/**
 * @file   material_non_local.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 28 11:17:41 2011
 *
 * @brief  Material class that handle the non locality of a law for example damage.
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
#include "aka_common.hh"
#include "material.hh"
#include "aka_grid.hh"
#include "fem.hh"

#include "weight_function.hh"

namespace akantu {
  class GridSynchronizer;
}

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_NON_LOCAL_HH__
#define __AKANTU_MATERIAL_NON_LOCAL_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt dim, template <UInt> class WeightFunction = BaseWeightFunction>
class MaterialNonLocal : public virtual Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialNonLocal(SolidMechanicsModel & model, const ID & id = "");
  virtual ~MaterialNonLocal();

  template<typename T>
  class PairList : public ByElementType<ByElementTypeVector<T> > {};
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// read properties
  virtual bool parseParam(const std::string & key, const std::string & value,
			  const ID & id);

  /// initialize the material computed parameter
  virtual void initMaterial();

  // void computeQuadraturePointsNeighborhoudVolumes(ByElementTypeReal & volumes) const;

  virtual void updateResidual(GhostType ghost_type);

  virtual void computeAllNonLocalStresses(GhostType ghost_type = _not_ghost);

  // void removeDamaged(const ByElementTypeReal & damage, Real thresold);

  void savePairs(const std::string & filename) const;
  void neighbourhoodStatistics(const std::string & filename) const;

protected:
  void updatePairList(const ByElementTypeReal & quadrature_points_coordinates);

  void computeWeights(const ByElementTypeReal & quadrature_points_coordinates);

  void createCellList(ByElementTypeReal & quadrature_points_coordinates);

  void fillCellList(const ByElementTypeReal & quadrature_points_coordinates,
		    const GhostType & ghost_type);

  /// constitutive law
  virtual void computeNonLocalStresses(GhostType ghost_type = _not_ghost) = 0;

  template<typename T>
  void weightedAvergageOnNeighbours(const ByElementTypeVector<T> & to_accumulate,
				    ByElementTypeVector<T> & accumulated,
				    UInt nb_degree_of_freedom,
				    GhostType ghost_type2 = _not_ghost) const;


  // template<typename T>
  // void accumulateOnNeighbours(const ByElementTypeVector<T> & to_accumulate,
  // 			      ByElementTypeVector<T> & accumulated,
  // 			      UInt nb_degree_of_freedom) const;


  virtual inline UInt getNbDataToPack(const Element & element,
				      SynchronizationTag tag) const;

  virtual inline UInt getNbDataToUnpack(const Element & element,
					SynchronizationTag tag) const;

  virtual inline void packData(CommunicationBuffer & buffer,
			       const Element & element,
			       SynchronizationTag tag) const;

  virtual inline void unpackData(CommunicationBuffer & buffer,
				 const Element & element,
				 SynchronizationTag tag);

  virtual UInt getNbDataToPack(SynchronizationTag tag) const { return Material::getNbDataToPack(tag); }
  virtual UInt getNbDataToUnpack(SynchronizationTag tag) const { return Material::getNbDataToUnpack(tag); }

  virtual void packData(CommunicationBuffer & buffer,
			const UInt index,
			SynchronizationTag tag) const {
    Material::packData(buffer, index, tag);
  }

  virtual void unpackData(CommunicationBuffer & buffer,
			const UInt index,
			SynchronizationTag tag) {
    Material::unpackData(buffer, index, tag);
  }


  virtual inline void onElementsAdded(const Vector<Element> & element_list);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  void registerNonLocalVariable(ByElementTypeReal & local,
				ByElementTypeReal & non_local,
				UInt nb_degree_of_freedom) {
    ID id = local.getID();
    NonLocalVariable & non_local_variable = non_local_variables[id];

    non_local_variable.local_variable = &local;
    non_local_variable.non_local_variable = &non_local;
    non_local_variable.non_local_variable_nb_component = nb_degree_of_freedom;
  }

  AKANTU_GET_MACRO(PairList, pair_list, const PairList<UInt> &)

  AKANTU_GET_MACRO(Radius, radius, Real);

  AKANTU_GET_MACRO(CellList, *cell_list, const RegularGrid<QuadraturePoint> &)

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the non local radius
  Real radius;

  /// the weight function used
  WeightFunction<dim> * weight_func;

private:
  /// the pairs of quadrature points
  PairList<UInt> pair_list;
  /// the weights associated to the pairs
  PairList<Real> pair_weight;

  /// the regular grid to construct/update the pair lists
  RegularGrid<QuadraturePoint> * cell_list;

  /// the types of the existing pairs
  typedef std::set< std::pair<ElementType, ElementType> > pair_type;
  pair_type existing_pairs[2];

  /// specify if the weights should be updated and at which rate
  UInt update_weights;

  /// count the number of calls of computeStress
  UInt compute_stress_calls;

  struct NonLocalVariable {
    ByElementTypeVector<Real> * local_variable;
    ByElementTypeVector<Real> * non_local_variable;
    UInt non_local_variable_nb_component;
  };

  std::map<ID, NonLocalVariable> non_local_variables;

  bool is_creating_grid;

  GridSynchronizer * grid_synchronizer;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_non_local_inline_impl.cc"


__END_AKANTU__

#endif /* __AKANTU_MATERIAL_NON_LOCAL_HH__ */
