/**
 * @file   non_local_neighborhood_base.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Sep 21 15:43:26 2015
 *
 * @brief  Non-local neighborhood base class
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
#ifndef __AKANTU_NON_LOCAL_NEIGHBORHOOD_BASE_HH__
#define __AKANTU_NON_LOCAL_NEIGHBORHOOD_BASE_HH__
/* -------------------------------------------------------------------------- */
#include "neighborhood_base.hh"
#include "parsable.hh"
/* -------------------------------------------------------------------------- */



__BEGIN_AKANTU__

class NonLocalNeighborhoodBase : public NeighborhoodBase,
				 public Parsable{
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  NonLocalNeighborhoodBase(const SolidMechanicsModel & model, 
			   const ElementTypeMapReal & quad_coordinates,
			   const ID & id = "non_local_neighborhood",
			   const MemoryID & memory_id = 0);
  virtual ~NonLocalNeighborhoodBase();


  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:

  /// create grid synchronizer and exchange ghost cells
  virtual void createGridSynchronizer();

  /// compute weights, for instance needed for non-local damage computation
  virtual void computeWeights() {};

  /// compute the non-local counter part for a given element type map
  virtual void weightedAverageOnNeighbours(const ElementTypeMapReal & to_accumulate,
					   ElementTypeMapReal & accumulated,
					   UInt nb_degree_of_freedom,
					   const GhostType & ghost_type2) const {};

  /// update the weights for the non-local averaging
  virtual void updateWeights() {};

  /// register a new non-local variable in the neighborhood
  virtual void registerNonLocalVariable(const ID & id) {};

  /// clean up the unneccessary ghosts
  void cleanupExtraGhostElements(std::set<Element> & relevant_ghost_elements);


protected:

  /// create the grid
  void createGrid();

/* -------------------------------------------------------------------------- */
/* DataAccessor inherited members                                             */
/* -------------------------------------------------------------------------- */
public:

  virtual inline UInt getNbDataForElements(const Array<Element> & elements,
					     SynchronizationTag tag) const {return 0; };

  virtual inline void packElementData(CommunicationBuffer & buffer,
  				      const Array<Element> & elements,
  				      SynchronizationTag tag) const {};

  virtual inline void unpackElementData(CommunicationBuffer & buffer,
  					const Array<Element> & elements,
  					SynchronizationTag tag) {};

/* -------------------------------------------------------------------------- */
/* Accessors                                                                  */
/* -------------------------------------------------------------------------- */
public:
  AKANTU_GET_MACRO(NonLocalVariables, non_local_variables, const std::set<ID> &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// list of non-local variables associated to the neighborhood
  std::set<ID> non_local_variables;
};



__END_AKANTU__
#endif /* __AKANTU_NON_LOCAL_NEIGHBORHOOD_BASE_HH__ */
