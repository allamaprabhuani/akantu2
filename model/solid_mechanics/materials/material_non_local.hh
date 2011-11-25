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

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_NON_LOCAL_HH__
#define __AKANTU_MATERIAL_NON_LOCAL_HH__

__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
class BaseWeightFonction {
public:
  BaseWeightFonction(Real radius,
                 __attribute__((unused)) ElementType type1,
                 __attribute__((unused)) GhostType ghost_type1,
                 __attribute__((unused)) ElementType type2,
                 __attribute__((unused)) GhostType ghost_type2) :
    radius(radius), r2(radius*radius) {
  }

  /* ------------------------------------------------------------------------ */
  inline Real operator()(Real r,
                         __attribute__((unused)) UInt q1,
                         __attribute__((unused)) UInt q2) {
    Real weight = 0;
    if(r <= radius) {
      Real alpha = (1. - r*r * r2);
      weight = alpha * alpha;
      //	*weight = 1 - sqrt(r / radius);
    }
    return weight;
  }

protected:
  Real radius;
  Real r2;
};



/* -------------------------------------------------------------------------- */


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

  template<class WeightFonction>
  void computeWeights(const ByElementTypeReal & quadrature_points_coordinates);

  void computeQuadraturePointsNeighborhoudVolumes(ByElementTypeReal & volumes) const;

  template<typename T>
  void accumulateOnNeighbours(const ByElementTypeVector<T> & to_accumulate,
			      ByElementTypeVector<T> & accumulated,
			      UInt nb_degree_of_freedom) const;

  template<typename T>
  void weigthedAvergageOnNeighbours(const ByElementTypeVector<T> & to_accumulate,
				    ByElementTypeVector<T> & accumulated,
				    UInt nb_degree_of_freedom) const;

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  // Real getStableTimeStep(Real h, const Element & element = ElementNull) {
  //   AKANTU_DEBUG_TO_IMPLEMENT();
  //   return 0.;
  // };

  // void computeStress(ElementType el_type,
  // 		     GhostType ghost_type = _not_ghost) {
  //   AKANTU_DEBUG_TO_IMPLEMENT();
  // };

  virtual void updateResidual(Vector<Real> & displacement, GhostType ghost_type);

  /// constitutive law
  virtual void computeNonLocalStress(Vector<Real> & averaged_field,
				     ElementType el_type,
				     GhostType ghost_type = _not_ghost) = 0;

  virtual void computeNonLocalStress(GhostType ghost_type = _not_ghost) = 0;

  void removeDamaged(const ByElementTypeReal & damage, Real thresold);

  // void computeTangentStiffness(const ElementType & el_type,
  // 			       Vector<Real> & tangent_matrix,
  // 			       GhostType ghost_type = _not_ghost) {
  //   AKANTU_DEBUG_TO_IMPLEMENT();
  // };


  void savePairs(const std::string & filename) const;


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(PairList, pair_list, const PairList<UInt> &)

  AKANTU_GET_MACRO(Radius, radius, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  PairList<UInt> pair_list;
  PairList<Real> pair_weigth;

  Real radius;

  //  ByElementTypeReal quadrature_points_coordinates;

  std::set< std::pair<ElementType, ElementType> > existing_pairs;
};


/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "material_non_local_inline_impl.cc"
#endif

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MaterialNonLocal & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_MATERIAL_NON_LOCAL_HH__ */
