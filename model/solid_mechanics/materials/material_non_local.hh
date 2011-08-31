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

class MaterialNonLocal : public Material {
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

  void updatePairList();

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {};

  virtual void computeWeights();

  void computeQuadraturePointsNeighborhoudVolumes(ByElementTypeReal & volumes) const;

  template<typename T>
  void accumulateOnNeighbours(const ByElementTypeVector<T> & to_accumulate,
			      ByElementTypeVector<T> & accumulated,
			      UInt nb_degree_of_freedom) const;

  template<typename T>
  void weigthedAvergageOnNeighbours(const ByElementTypeVector<T> & to_accumulate,
				    ByElementTypeVector<T> & accumulated,
				    UInt nb_degree_of_freedom) const;

  Real getStableTimeStep(Real h, const Element & element = ElementNull) { return 0.; };

  void computeStress(ElementType el_type,
		     GhostType ghost_type = _not_ghost) {};

  void computeTangentStiffness(const ElementType & el_type,
			       Vector<Real> & tangent_matrix,
			       GhostType ghost_type = _not_ghost) {};


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

  ByElementTypeReal quadrature_points_coordinates;

  std::set< std::pair<ElementType, ElementType> > existing_pairs;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_non_local_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const MaterialNonLocal & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_MATERIAL_NON_LOCAL_HH__ */
