/**
 * @file   fragment_manager.hh
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Mon Jan 20 14:26:42 2014
 *
 * @brief  Group manager to handle fragments
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

#ifndef __AKANTU_FRAGMENT_MANAGER_HH__
#define __AKANTU_FRAGMENT_MANAGER_HH__

#include "group_manager.hh"
#include "solid_mechanics_model_cohesive.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */

class FragmentManager : public GroupManager {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  FragmentManager(SolidMechanicsModelCohesive & model,
		  bool dump_data = true,
		  const ID & id = "fragment_manager",
		  const MemoryID & memory_id = 0);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
private:

  /// store mass density per quadrature point
  void storeMassDensityPerQuadraturePoint();

  /// integrate an elemental field multiplied by density on global
  /// fragments
  void integrateFieldOnFragments(ByElementTypeReal & field,
				 Array<Real> & output);

  /// compute fragments' mass
  void computeMass();

  /// create dump data for a single array
  template <typename T>
  void createDumpDataArray(Array<T> & data, std::string name);

public:

  /// build fragment list
  void buildFragments();

  /// compute fragments' center of mass
  void computeCenterOfMass();

  /// compute fragments' velocity
  void computeVelocity();

  /// computes principal moments of inertia with respect to the center
  /// of mass of each fragment
  void computeInertiaMoments();

  /// compute all fragments' data
  void computeAllData();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// get number of fragments
  AKANTU_GET_MACRO(NbFragment, global_nb_fragment, UInt);

  /// get fragments' mass
  AKANTU_GET_MACRO(Mass, mass, const Array<Real> &);

  /// get fragments' center of mass
  AKANTU_GET_MACRO(CenterOfMass, mass_center, const Array<Real> &);

  /// get fragments' velocity
  AKANTU_GET_MACRO(Velocity, velocity, const Array<Real> &);

  /// get fragments' principal moments of inertia
  AKANTU_GET_MACRO(MomentsOfInertia, inertia_moments, const Array<Real> &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// local_fragment index
  Array<UInt> fragment_index;

  /// global number of fragments (parallel simulations)
  UInt global_nb_fragment;

  /// number of fragments
  UInt nb_fragment;

  /// cohesive solid mechanics model associated
  SolidMechanicsModelCohesive & model;

  /// fragments' center of mass
  Array<Real> mass_center;

  /// fragments' mass
  Array<Real> mass;

  /// fragments' velocity
  Array<Real> velocity;

  /// fragments' principal moments of inertia with respect to the
  /// center of mass
  Array<Real> inertia_moments;

  /// quadrature points' coordinates
  ByElementTypeReal quad_coordinates;

  /// mass density per quadrature point
  ByElementTypeReal mass_density;

  /// dump data
  bool dump_data;

};

__END_AKANTU__

#endif /* __AKANTU_FRAGMENT_MANAGER_HH__ */
