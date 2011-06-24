/**
 * @file   heat_transfer_model.cc
 * @author Rui WANG<rui.wang@epfl.ch>
 * @date   Fri May  4 13:46:43 2011
 *
 * @brief  Implementation of HeatTransferModel class
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
 *el_size
 */

/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
#include "heat_transfer_model.hh"
#include "aka_math.hh"
//#include "integration_scheme_2nd_order.hh"
#include "aka_common.hh"
#include <fstream>
#include <iostream>
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "static_communicator.hh"
#include "sparse_matrix.hh"
#include "solver.hh"
//#include "solid_mechanics_model.hh"
#include"sparse_matrix.hh"
// #include "fem_template.cc"
 #include "fem_template.hh"
using namespace std;
#ifdef AKANTU_USE_MUMPS
#include "solver_mumps.hh"
#endif
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

HeatTransferModel::HeatTransferModel(Mesh & mesh,
 		      UInt dim,
 		      const ModelID & id,
		      const MemoryID & memory_id):
  Model(id, memory_id),
  spatial_dimension(dim) 
{
  AKANTU_DEBUG_IN();

  createSynchronizerRegistry(this);

  if (spatial_dimension == 0) spatial_dimension = mesh.getSpatialDimension();
  registerFEMObject<MyFEMType>("HeatTransferFEM",mesh,spatial_dimension);

 this->temperature= NULL;
 for (UInt t = _not_defined; t < _max_element_type; ++t) {
   this->temperature_gradient[t] = NULL;   
 }

 this->heat_flux = NULL;
 this->boundary = NULL;

 AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* Initialisation                                                             */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initModel()
{
  getFEM().initShapeFunctions(_not_ghost);
  getFEM().initShapeFunctions(_ghost);
}
/* -------------------------------------------------------------------------- */
void HeatTransferModel::initPBC(UInt x, UInt y, UInt z){
  Model::initPBC(x,y,z);
  PBCSynchronizer * synch = new PBCSynchronizer(pbc_pair);
  synch_registry->registerSynchronizer(*synch,akantu::_gst_htm_capacity);
  synch_registry->registerSynchronizer(*synch,akantu::_gst_htm_temperature);
  changeLocalEquationNumberforPBC(pbc_pair,1);
}
/* -------------------------------------------------------------------------- */
void HeatTransferModel::initVectors()
{
  AKANTU_DEBUG_IN();

  UInt nb_nodes = getFEM().getMesh().getNbNodes();
  std::stringstream sstr_temp; sstr_temp << id << ":temperature";
  std::stringstream sstr_heat_flux; sstr_heat_flux << id << ":heat flux";
  std::stringstream sstr_lump; sstr_lump<<id<<":lumped";
  std::stringstream sstr_boun; sstr_boun<<id<<":boundary";
  temperature = &(alloc<Real>(sstr_temp.str(), nb_nodes, 1, REAL_INIT_VALUE));
  heat_flux= &(alloc<Real>(sstr_heat_flux.str(), nb_nodes, 1, REAL_INIT_VALUE));
  capacity_lumped= &(alloc<Real>(sstr_lump.str(), nb_nodes, 1, REAL_INIT_VALUE));
  boundary= &(alloc<bool>(sstr_boun.str(), nb_nodes, 1, false));
  const Mesh::ConnectivityTypeList & type_list = 
    this->getFEM().getMesh().getConnectivityTypeList(_not_ghost);
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it)
    {
      if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;
      std::stringstream sstr_temp_grad; 
      sstr_temp_grad << id << ":temperature gradient:" << *it;
      UInt nb_element = getFEM().getMesh().getNbElement(*it);
      UInt nb_quad_points = this->getFEM().getNbQuadraturePoints(*it) * nb_element;
      std::cout << "Allocating " << sstr_temp_grad.str() << std::endl;
      temperature_gradient[*it] = &(alloc<Real>(sstr_temp_grad.str(), 
						nb_quad_points, 
						spatial_dimension, REAL_INIT_VALUE)); 
    }

  dof_synchronizer = new DOFSynchronizer(getFEM().getMesh(),1);
  dof_synchronizer->initLocalDOFEquationNumbers();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

HeatTransferModel::~HeatTransferModel() 
{
  AKANTU_DEBUG_IN();
 
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

//set boundary condition
void setBoundaryCondition()
{
  AKANTU_DEBUG_IN();
  // Mesh::getcomputeBoundingBox();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void HeatTransferModel::assembleCapacityLumped(const ElementType &el_type)
{
   AKANTU_DEBUG_IN();
   
   FEM & fem = getFEM();
  
  UInt nb_quadrature_points = getFEM().getNbQuadraturePoints(el_type);
  UInt nb_element = getFEM().getMesh().getNbElement(el_type);
  //???
  Vector<Real> rho_1 (nb_element * nb_quadrature_points,1, capacity * density);
 
  const Mesh::ConnectivityTypeList & type_list = fem.getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it)
  {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;
    
    fem.assembleFieldLumped(rho_1,1,*capacity_lumped,
			    dof_synchronizer->getLocalDOFEquationNumbers(),
			    *it, _not_ghost);
  }
  
  getSynchronizerRegistry().synchronize(akantu::_gst_htm_capacity);

  AKANTU_DEBUG_OUT();
}
  
/* -------------------------------------------------------------------------- */


// compute the heat flux 
 void HeatTransferModel::updateHeatFlux()
 {
   AKANTU_DEBUG_IN();
   heat_flux->clear();

   const Mesh::ConnectivityTypeList & type_list = this->getFEM().getMesh().getConnectivityTypeList(_not_ghost);
   Mesh::ConnectivityTypeList::const_iterator it;
 
   for(it = type_list.begin(); it != type_list.end(); ++it)
   {
     if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

     const Vector<Real> * shapes_derivatives;
      shapes_derivatives = &(this->getFEM().getShapesDerivatives(*it));

     UInt nb_element = getFEM().getMesh().getNbElement(*it);     
     UInt nb_quadrature_points = getFEM().getNbQuadraturePoints(*it);
     UInt size_of_shapes_derivatives = shapes_derivatives->getNbComponent();
     UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);

     /*problem is included in computeTemperatureGradient*/
     //const Vector<Real> *gT = computeTemperatureGradient(*it);
     computeTemperatureGradient(*it);

     Real * gT_val = temperature_gradient[*it]->values;
     Real * k_gT_val = new Real[spatial_dimension];
     Real * shapes_derivatives_val = shapes_derivatives->values;
     UInt bt_k_gT_size = nb_nodes_per_element;
     Vector<Real> * bt_k_gT = new Vector <Real> (nb_quadrature_points*nb_element, bt_k_gT_size);
     Real * bt_k_gT_val = bt_k_gT->values;
     
     for (UInt el = 0; el < nb_element; ++el)
     {
       for (UInt i = 0; i < nb_quadrature_points; ++i)
	 { 
	   Math::matrixt_vector(spatial_dimension, spatial_dimension,
				conductivity,
				gT_val,
				k_gT_val);
	   gT_val += spatial_dimension;

	   Math::matrix_vector(nb_nodes_per_element, spatial_dimension,
			       shapes_derivatives_val,
			       k_gT_val,
			       bt_k_gT_val);
	 
	   shapes_derivatives_val += nb_nodes_per_element * spatial_dimension;
	   bt_k_gT_val += bt_k_gT_size;
	 }
     }     
     
     Vector<Real> * q_e = new Vector<Real>(0,bt_k_gT_size,"q_e");
     this->getFEM().integrate(*bt_k_gT, *q_e, bt_k_gT_size, *it);
     delete bt_k_gT;
     
     // SparseMatrix & K = new SparseMatrix(this->getFEM().getMesh().getNbNodes());
     // const Vector<Int> & equation_number = this->getEquationNumber();
     // this->getFEM().assembleMatrix(*q_e, K, equation_number, spatial_dimension, *it);

     //the previous version of the assemble Vector, it is right, just ignore the external heat flux
     this->getFEM().assembleVector(*q_e, *heat_flux, 
				   dof_synchronizer->getLocalDOFEquationNumbers(),
				   1, *it,_not_ghost,NULL,-1);
     delete q_e;
     
     AKANTU_DEBUG_OUT();
   }
}
 

/* -------------------------------------------------------------------------- */



 //compute temperature gradient

void HeatTransferModel:: computeTemperatureGradient(const ElementType &el_type )
{
  AKANTU_DEBUG_IN();
 this->getFEM().gradientOnQuadraturePoints(*temperature, *temperature_gradient[el_type], 1 ,el_type);



  AKANTU_DEBUG_OUT();
}

/*except updateflux, everything is fine*/
//compute the temperature 
void HeatTransferModel::updateTemperature()
{
  AKANTU_DEBUG_IN();
  UInt nb_nodes=getFEM().getMesh().getNbNodes();

  for(UInt i=0; i < nb_nodes; i++)
  {
    if(!((*boundary)(i)))
    {	
      (*temperature)(i,0) += (*heat_flux)(i,0) * time_step / (*capacity_lumped)(i,0);
      (*temperature)(i,0) = std::max((*temperature)(i,0), 0.0);
    }      
  }

  synch_registry->synchronize(akantu::_gst_htm_temperature);
  AKANTU_DEBUG_OUT();
}



/* -------------------------------------------------------------------------- */


//choose scheme for temperature iteration
/*a=0.0, the forward Euler; a=0.5, the Crank-Nicolson;a=2.0/3, Galerkin;a=1.0,backward Euler*/
void HeatTransferModel::integrationScheme1stOrder(Real thelta, UInt N,  Vector<Real> * temperature)
{
  AKANTU_DEBUG_IN();
  
 //remember to use time_step
  if(thelta==0.0)
  {
      // temperature->values[i] +=  temperature->values[i] + heat_flux->values[i]/lumped->values[i];
  }
  else if(thelta==0.5)
  {

  }
 
  else if(thelta==1.0)
  {
      
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getStableTimeStep()
{
  AKANTU_DEBUG_IN();
 
  Real conductivitymax = -std::numeric_limits<Real>::max();

  Real el_size;
  
  Real min_el_size = std::numeric_limits<Real>::max();

  Real * coord    = getFEM().getMesh().getNodes().values;
 

  const Mesh::ConnectivityTypeList & type_list = getFEM().getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it)
 {
    if(getFEM().getMesh().getSpatialDimension(*it) != spatial_dimension) continue;
    
    UInt nb_nodes_per_element = getFEM().getMesh().getNbNodesPerElement(*it);
    UInt nb_element           = getFEM().getMesh().getNbElement(*it);

    UInt * conn         = getFEM().getMesh().getConnectivity(*it).values;
    Real * u = new Real[nb_nodes_per_element*spatial_dimension];

    for (UInt el = 0; el < nb_element; ++el) 
    {
      UInt el_offset  = el * nb_nodes_per_element;
      for (UInt n = 0; n < nb_nodes_per_element; ++n) 
      {
	UInt offset_conn = conn[el_offset + n] * spatial_dimension;
	memcpy(u + n * spatial_dimension,
	       coord + offset_conn,
	       spatial_dimension * sizeof(Real));


      }

     el_size    = getFEM().getElementInradius(u, *it);
     
     //get the biggest parameter from k11 until k33//
     for(UInt i = 0; i < spatial_dimension * spatial_dimension; i++)
     {
       if(conductivity[i] > conductivitymax)
	 conductivitymax=conductivity[i];
     }  
     min_el_size = std::min(min_el_size, el_size);
     

    }
    cout<<min_el_size<<endl;
    cout<<"the minimum size:"<<conductivitymax<<endl;
    delete [] u;
  }


  Real min_dt    = 2* min_el_size * min_el_size * density * capacity/conductivitymax ;
  time_step = min_dt;
  return time_step;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void HeatTransferModel::readMaterials(const std::string & filename) {
  Parser parser;
  parser.open(filename);
  std::string mat_type = parser.getNextSection("heat");
  UInt mat_count = 0;

  if (mat_type != ""){
    parser.readSection<HeatTransferModel>(*this);
  }
  else
    AKANTU_DEBUG_ERROR("did not find any section with material info for heat conduction");
}

/* -------------------------------------------------------------------------- */

void HeatTransferModel::setParam(const std::string & key, const std::string & value) {
  std::stringstream str(value);
  if (key == "conductivity"){
    conductivity = new Real[spatial_dimension * spatial_dimension];
    for(int i=0;i<3;i++)
      for(int j=0; j<3;j++)
	{
	  if (i< spatial_dimension && j < spatial_dimension){
	    str >> conductivity[i*spatial_dimension+j];
	    AKANTU_DEBUG_INFO("conductivity(" << i << "," << j << ") = " 
			      << conductivity[i*spatial_dimension+j]);
	  }
	  else {
	    Real tmp;
	    str >> tmp;
	  }
	}
  }
  else if (key == "capacity"){
    str >> capacity;
    AKANTU_DEBUG_INFO("The capacity of the material is:" << capacity);
  }
  else if (key == "density"){
    str >> density;    
    AKANTU_DEBUG_INFO("The density of the material is:" << density);
  }	
}
/* -------------------------------------------------------------------------- */


__END_AKANTU__
