/**
 * @file   test_cohesive.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Fri Feb 24 14:32:31 2012
 *
 * @brief  Test for cohesive elements
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

#include <limits>
#include <fstream>
#include <iostream>


/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "material.hh"
#if defined(AKANTU_USE_IOHELPER)
#  include "io_helper.hh"
#endif
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  initialize(argc, argv);

  debug::setDebugLevel(dblError);

  const UInt spatial_dimension = 2;
  const ElementType type = _triangle_6;

  Mesh mesh(spatial_dimension);
  MeshIOMSH mesh_io;
  mesh_io.read("simple.msh", mesh);

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull("material.dat", _static);
  
  const Mesh & mesh_facets = model.getMeshFacets();

  const ElementType type_facet = mesh.getFacetElementType(type);
  UInt nb_facet = mesh_facets.getNbElement(type_facet);
  Vector<UInt> facet_insertion;
  Real * bary_facet = new Real[spatial_dimension];
  for (UInt f = 0; f < nb_facet; ++f) {
    mesh_facets.getBarycenter(f, type_facet, bary_facet);
    if (bary_facet[1] == 1) facet_insertion.push_back(f);
  }
  delete[] bary_facet;

  model.insertCohesiveElements(facet_insertion);

  /// boundary conditions
  Vector<bool> & boundary = model.getBoundary();
  UInt nb_nodes = mesh.getNbNodes();
  UInt nb_element = mesh.getNbElement(type);
  Vector<Real> &position = const_cast <Vector<Real>&> (mesh.getNodes());
  Vector<Real> & displacement = model.getDisplacement();
  
  
  for (UInt n = 0; n < nb_nodes; ++n) {
    
    if (std::abs(position(n,1))< Math::getTolerance()){
      boundary(n, 1) = true;
      displacement(n,1) = 0.0;
    }

    if ((std::abs(position(n,0))< Math::getTolerance())&& (position(n,1)< 1.1)){
      boundary(n, 0) = true;
      displacement(n,0) = 0.0;
    }
    if ((std::abs(position(n,0)-1)< Math::getTolerance())&&(std::abs(position(n,1)-1)< Math::getTolerance())){
      boundary(n, 0) = true;
      displacement(n,0) = 0.0;
    }

    if (std::abs(position(n,1)-2)< Math::getTolerance()){
      boundary(n, 1) = true;
    }
  }
  
#if defined(AKANTU_USE_IOHELPER)  
  iohelper::ElemType paraview_type = iohelper::TRIANGLE2;

  /// initialize the paraview output
  iohelper::DumperParaview dumper;
  dumper.SetPoints(model.getFEM().getMesh().getNodes().values,
  		   spatial_dimension, nb_nodes, "implicit");
  dumper.SetConnectivity((int *)model.getFEM().getMesh().getConnectivity(type).values,
  			 paraview_type, nb_element, iohelper::C_MODE);
  dumper.AddNodeDataField(model.getDisplacement().values,
  			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model.getVelocity().values,
  			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model.getAcceleration().values,
  			  spatial_dimension, "acceleration");
  dumper.AddNodeDataField(model.getForce().values,
  			  spatial_dimension, "applied_force");
  dumper.AddNodeDataField(model.getResidual().values,
   			  spatial_dimension, "forces");
  dumper.AddElemDataField(model.getMaterial(0).getStrain(type).values,
  			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(model.getMaterial(0).getStress(type).values,
  			  spatial_dimension*spatial_dimension, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("applied_force", 1);
  dumper.SetEmbeddedValue("forces", 1);
  dumper.SetPrefix("paraview");
  dumper.Init();
  dumper.Dump();
#endif

  const MaterialCohesive & mat_coh = dynamic_cast< const MaterialCohesive &> (model.getMaterial(1));

  const Vector<Real> & opening = mat_coh.getOpening(_cohesive_2d_6);
  //const Vector<Real> & traction = mat_coh.getTraction(_cohesive_2d_6);

  //model.getStiffnessMatrix().saveMatrix("K.mtx");
  model.updateResidual();
  const Vector<Real> & residual = model.getResidual();

  UInt max_step = 1000;
  Real increment =  3./max_step;
  Real error_tol = 10e-6;

  std::ofstream fout;
  fout.open("output");

  /// Main loop
  for ( UInt nstep = 0; nstep < max_step; ++nstep){ 
    Real norm = 10;
    UInt count = 0;
    
    for (UInt n = 0; n < nb_nodes; ++n) {
      if (std::abs(position(n,1)-2)< Math::getTolerance()){
	displacement(n,1) += increment;
      }
    }

    do{
      std::cout << "Iter : " << ++count << " - residual norm : " << norm << std::endl;
      model.assembleStiffnessMatrix();
      model.solveStatic();
      model.updateResidual();
     
    } while(!model.testConvergenceResidual(1e-5, norm) && (count < 100))  ;

    std::cout << "Step : " << nstep << " - residual norm : " << norm << std::endl;

#if defined(AKANTU_USE_IOHELPER)
    dumper.Dump(); 
#endif

    Real resid = 0;
    for (UInt n = 0; n < nb_nodes; ++n) {
      if (std::abs(position(n,1)-2)< Math::getTolerance()){
	resid += residual.values[spatial_dimension* n + 1];
      }
    }

    Real analytical = exp(1) * std::abs(opening.values[3]) * exp (-std::abs(opening.values[3])/0.5)/0.5;

    // the residual force is comparing with the theoretical value of the cohesive law
    error_tol  = std::abs((- resid - analytical)/analytical);

    fout << nstep << " " << -resid << " " << analytical << " " << error_tol << std::endl;
 
    if (error_tol > 2e-5) 
      return EXIT_FAILURE;
  }

  fout.close();

  finalize();

  return EXIT_SUCCESS;
}



