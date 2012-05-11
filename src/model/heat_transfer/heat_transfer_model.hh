/**
 * @file   heat_transfer_model.hh
 * @author Rui WANG<rui.wang@epfl.ch>
 * @date   Fri May  4 13:35:55 2011
 *
 * @brief  Model of Heat Transfer
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

#ifndef __AKANTU_HEAT_TRANSFER_MODEL_HH__
#define __AKANTU_HEAT_TRANSFER_MODEL_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_memory.hh"
#include "model.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"

namespace akantu {
  class IntegrationScheme1stOrder;
//    class Solver;
//    class SparseMatrix;
}

__BEGIN_AKANTU__

class HeatTransferModel : public Model, public DataAccessor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef FEMTemplate<IntegratorGauss,ShapeLagrange> MyFEMType;

  HeatTransferModel(UInt spatial_dimension,
		    const ID & id = "heat_transfer_model",
  		    const MemoryID & memory_id = 0) ;

  HeatTransferModel(Mesh & mesh,
		    UInt spatial_dimension = 0,
		    const ID & id = "heat_transfer_model",
		    const MemoryID & memory_id = 0);

  virtual ~HeatTransferModel() ;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:

  /// generic function to initialize everything ready for explicit dynamics
  void initFull(const std::string & material_file);

  /// set the parameters
  bool setParam(const std::string & key, const std::string & value);

  /// read one material file to instantiate all the materials
  void readMaterials(const std::string & filename);

  /// allocate all vectors
  void initVectors();

  /// register the tags associated with the parallel synchronizer
  void initParallel(MeshPartition * partition, DataAccessor * data_accessor=NULL);

  /// initialize the model
  void initModel();

  /// init PBC synchronizer
  void initPBC();

  /// function to print the contain of the class
  virtual void printself(__attribute__ ((unused)) std::ostream & stream,
			 __attribute__ ((unused)) int indent = 0) const {};

  /* ------------------------------------------------------------------------ */
  /* Methods for explicit                                                     */
  /* ------------------------------------------------------------------------ */
public:

  /// compute and get the stable time step
  Real getStableTimeStep();

  /// compute the heat flux
  void updateResidual();

  /// calculate the lumped capacity vector for heat transfer problem
  void assembleCapacityLumped();

  /// update the temperature from the temperature rate
  void explicitPred();

  /// update the temperature rate from the increment
  void explicitCorr();


  // /// initialize the heat flux
  // void initializeResidual(Vector<Real> &temp);
  // /// initialize temperature
  // void initializeTemperature(Vector<Real> &temp);

private:

  /// solve the system in temperature rate  @f$C\delta \dot T = q_{n+1} - C \dot T_{n}@f$
  void solveExplicitLumped();

  /// compute the heat flux on ghost types
  void updateResidual(const GhostType & ghost_type);

  /// calculate the lumped capacity vector for heat transfer problem (w ghosttype)
  void assembleCapacityLumped(const GhostType & ghost_type);
  
  /// compute the conductivity tensor for each quadrature point in an array
  void computeConductivityOnQuadPoints(const GhostType & ghost_type);

  /// compute vector k \grad T for each quadrature point 
  void computeKgradT(const GhostType & ghost_type);


  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:

  inline UInt getNbDataToPack(const Element & element,
			      SynchronizationTag tag) const;
  inline UInt getNbDataToUnpack(const Element & element,
				SynchronizationTag tag) const;
  inline void packData(CommunicationBuffer & buffer,
		       const Element & element,
		       SynchronizationTag tag) const;
  inline void unpackData(CommunicationBuffer & buffer,
			 const Element & element,
			 SynchronizationTag tag);

  inline UInt getNbDataToPack(SynchronizationTag tag) const;
  inline UInt getNbDataToUnpack(SynchronizationTag tag) const;
  inline void packData(CommunicationBuffer & buffer,
		       const UInt index,
		       SynchronizationTag tag) const;
  inline void unpackData(CommunicationBuffer & buffer,
			 const UInt index,
			 SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// get the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);
  /// get the current value of the time step
  AKANTU_GET_MACRO(TimeStep, time_step, Real);
  /// set the value of the time step
  AKANTU_SET_MACRO(TimeStep, time_step, Real);
  /// get the assembled heat flux
  AKANTU_GET_MACRO(Residual, *residual, Vector<Real>&);
  /// get the lumped capacity
  AKANTU_GET_MACRO(CapacityLumped, * capacity_lumped, Vector<Real>&);
  /// get the boundary vector
  AKANTU_GET_MACRO(Boundary, * boundary, Vector<bool>&);
  /// get the temperature gradient
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(TemperatureGradient, temperature_gradient, Real);
  /// get the conductivity on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ConductivityOnQpoints, conductivity_on_qpoints, Real);
  /// get the conductivity on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(TemperatureOnQpoints, temperature_on_qpoints, Real);
  /// get the temperature
  AKANTU_GET_MACRO(Temperature, *temperature, Vector<Real> &);
  /// get the temperature derivative
  AKANTU_GET_MACRO(TemperatureRate, *temperature_rate, Vector<Real> &);
  /// get the equation number Vector<Int>
  AKANTU_GET_MACRO(EquationNumber, *equation_number, const Vector<Int> &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  IntegrationScheme1stOrder * integrator;

  /// time step
  Real time_step;

  /// temperatures array
  Vector<Real> * temperature;

  /// temperatures derivatives array
  Vector<Real> * temperature_rate;

  /// increment array (@f$\delta \dot T@f$ or @f$\delta T@f$)
  Vector<Real> * increment;

  /// the spatial dimension
  UInt spatial_dimension;

  /// the density
  Real density;

  /// the speed of the changing temperature
  ByElementTypeReal temperature_gradient;

  /// temperature field on quadrature points
  ByElementTypeReal temperature_on_qpoints;

  /// conductivity tensor on quadrature points
  ByElementTypeReal conductivity_on_qpoints;

  /// vector k \grad T on quad points
  ByElementTypeReal k_gradt_on_qpoints;

  /// vector \int \grad N k \grad T
  ByElementTypeReal int_bt_k_gT;

  /// vector \grad N k \grad T
  ByElementTypeReal bt_k_gT;

  //external flux vector
  Vector<Real> * external_flux;

  /// residuals array
  Vector<Real> * residual;

  /// position of a dof in the K matrix
  Vector<Int> * equation_number;

  //lumped vector
  Vector<Real> * capacity_lumped;

  /// boundary vector
  Vector<bool> * boundary;

  //realtime
  Real time;

  ///capacity
  Real capacity;

  //conductivity matrix
  Real* conductivity;

  //linear variation of the conductivity (for temperature dependent conductivity)
  Real conductivity_variation;

  // reference temperature for the interpretation of temperature variation
  Real t_ref;

  //the biggest parameter of conductivity matrix
  Real conductivitymax;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "heat_transfer_model_inline_impl.cc"
#endif

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const HeatTransferModel & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__



#endif /* __AKANTU_HEAT_TRANSFER_MODEL_HH__ */
