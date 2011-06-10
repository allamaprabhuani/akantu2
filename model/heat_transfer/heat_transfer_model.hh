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
#include <fstream>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "model.hh"
#include "material.hh"
#include "parser.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
#include "fem.hh"
#include "mesh.hh"

#include "aka_memory.hh"
#include "element_class.hh"
#include "sparse_matrix.hh"


// namespace akantu {
//   class Solver;
//   class SparseMatrix;
// }

__BEGIN_AKANTU__

class HeatTransferModel : public Model {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef FEMTemplate<IntegratorGauss,ShapeLagrange> MyFEMType;
  
  HeatTransferModel(UInt spatial_dimension,
		    const ModelID & id = "heat_transfer_model",
  		    const MemoryID & memory_id = 0) ;  

  HeatTransferModel(Mesh & mesh,
		    UInt spatial_dimension = 0,
		    const ModelID & id = "heat_transfer_model",
		    const MemoryID & memory_id = 0);

  virtual ~HeatTransferModel() ;
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:
  
  /// set the parameters 
  void setParam(const std::string & key, const std::string & value);


  /// read one material file to instantiate all the materials
  void readMaterials(const std::string & filename);

  /// allocate all vectors
  void initVectors();

  /// initialize the model
  void initModel();
  
  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {};
  
  /* ------------------------------------------------------------------------ */
  /* Ghost Synchronizer inherited members                                     */
  /* ------------------------------------------------------------------------ */
public:

  inline virtual UInt getNbDataToPack(const Element & element,
				      GhostSynchronizationTag tag) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };

  inline virtual UInt getNbDataToUnpack(const Element & element,
					GhostSynchronizationTag tag) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };

  inline virtual void packData(CommunicationBuffer & buffer,
			       const Element & element,
			       GhostSynchronizationTag tag) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };

  inline virtual void unpackData(CommunicationBuffer & buffer,
				 const Element & element,
				 GhostSynchronizationTag tag) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  };


 
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// get the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  //* all the implementation of time step ----------------------------------- */
  /// get the current value of the time step
  AKANTU_GET_MACRO(TimeStep, time_step, Real);
  /// set the value of the time step
  AKANTU_SET_MACRO(TimeStep, time_step, Real);
 /// set the value of the time step
  AKANTU_GET_MACRO(HeatFlux,* heat_flux, Vector<Real>&);
/// set the value of the time step
  AKANTU_GET_MACRO(Lumped, * lumped, Vector<Real>&);
  AKANTU_GET_MACRO(Boundary, * boundary, Vector<bool>&);
  /// get the SolidMechanicsModel::velocity vector
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(TemperatureGradient, temperature_gradient, Vector<Real> &);
 /// get the SolidMechanicsModel::velocity vector
  AKANTU_GET_MACRO(Temperature, *temperature, Vector<Real> &);
 
 /// get the equation number Vector<Int>
  AKANTU_GET_MACRO(EquationNumber, *equation_number, const Vector<Int> &);
  
 
  // AKANTU_GET_MACRO(HeatFlux, *heat_flux, Vector<Real> &);
/// compute the stable time step
  Real getStableTimeStep();



 //initialize the heat flux
  void initializeHeatFlux(Vector<Real> &temp);
  //initialize temperature
  void initializeTemperature(Vector<Real> &temp);
  //set boundary condition
  void setBoundaryCondition();
  //compute temperature gradient
  void computeTemperatureGradient(const ElementType &el_type);

   //compute the heat flux
  void updateHeatFlux();

  //compute the temperature 
  void updateTemperature();

  //put the scheme into iteration
  void integrationScheme1stOrder(Real thelta, UInt N, Vector<Real> * temperature);

  //calculate the capacity matrix of heat transfer problem
  void assembleCapacityLumped(const ElementType &el_type);

  





  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// time step
  Real time_step;

  /// temperatures array
  Vector<Real> * temperature ;

  /// the spatial dimension
  UInt spatial_dimension;

  /// the density
  Real density;

  /// the speed of the changing temperature
  ByElementTypeReal temperature_gradient;

  /// K*T internal flux vector
  Vector<Real> * heat_flux;

  //external flux vector

  Vector<Real> * externalFlux;
  /// residuals array
  Vector<Real> * residual;

  /// position of a dof in the K matrix
  Vector<Int> * equation_number;

  //lumped vector
  Vector<Real> * lumped;

  /// boundary vector
  Vector<bool> * boundary;

  //realtime
  Real time;
  ///capacity
  Real capacity;

  //conductivity matrix
  Real* conductivity;

  //the biggest parameter of conductivity matrix
  Real conductivitymax;

  
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "heat_transfer_model_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const HeatTransferModel & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__



#endif /* __AKANTU_HEAT_TRANSFER_MODEL_HH__ */
