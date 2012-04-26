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


/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_MATERIAL_NON_LOCAL_HH__
#define __AKANTU_MATERIAL_NON_LOCAL_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<class WeightFunction = BaseWeightFunction>
class MaterialNonLocal : public virtual Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialNonLocal(Model & model, const ID & id = "");
  virtual ~MaterialNonLocal();

  template<typename T>
  class PairList : public ByElementType<ByElementTypeVector<T> > {};
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// read properties
  virtual bool setParam(const std::string & key, const std::string & value,
			const ID & id);

  /// initialize the material computed parameter
  virtual void initMaterial();

  void updatePairList(const ByElementTypeReal & quadrature_points_coordinates);

  void computeWeights(const ByElementTypeReal & quadrature_points_coordinates);

  // void computeQuadraturePointsNeighborhoudVolumes(ByElementTypeReal & volumes) const;

  template<typename T>
  void weightedAvergageOnNeighbours(const ByElementTypeVector<T> & to_accumulate,
				    ByElementTypeVector<T> & accumulated,
				    UInt nb_degree_of_freedom) const;

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  virtual void updateResidual(Vector<Real> & displacement, GhostType ghost_type);

  /// constitutive law
  virtual void computeNonLocalStress(GhostType ghost_type = _not_ghost) = 0;

  // void removeDamaged(const ByElementTypeReal & damage, Real thresold);

  void savePairs(const std::string & filename) const;
  void neighbourhoodStatistics(const std::string & filename) const;

protected:
  void createCellList(const ByElementTypeReal & quadrature_points_coordinates);

  // template<typename T>
  // void accumulateOnNeighbours(const ByElementTypeVector<T> & to_accumulate,
  // 			      ByElementTypeVector<T> & accumulated,
  // 			      UInt nb_degree_of_freedom) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(PairList, pair_list, const PairList<UInt> &)

  AKANTU_GET_MACRO(Radius, radius, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the non local radius
  Real radius;

  /// the weight function used
  WeightFunction * weight_func;

private:
  /// the pairs of quadrature points
  PairList<UInt> pair_list;
  /// the weights associated to the pairs
  PairList<Real> pair_weight;

  /// the regular grid to construct/update the pair lists
  RegularGrid<QuadraturePoint> * cell_list;

  /// the types of the existing pairs
  std::set< std::pair<ElementType, ElementType> > existing_pairs;

  /// specify if the weights should be updated and at which rate
  UInt update_weigths;

  /// count the number of calls of computeStress
  UInt compute_stress_calls;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "material_non_local_inline_impl.cc"
#endif

/// standard output stream operator
template<class WeightFunction>
inline std::ostream & operator <<(std::ostream & stream, const MaterialNonLocal<WeightFunction> & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_MATERIAL_NON_LOCAL_HH__ */
