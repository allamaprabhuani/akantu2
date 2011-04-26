/**
 * @file   test_solid_mechanics_model.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 14:34:13 2010
 *
 * @brief  test of the class SolidMechanicsModel
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

//#define CHECK_STRESS
akantu::ElementType type = akantu::_triangle_6;
#ifdef AKANTU_USE_IOHELPER
  akantu::UInt paraview_type = TRIANGLE2;
#endif //AKANTU_USE_IOHELPER

akantu::SolidMechanicsModel * model;
akantu::UInt spatial_dimension = 2;
akantu::UInt nb_nodes;
akantu::UInt nb_element;
akantu::UInt nb_quadrature_points;

akantu::Vector<akantu::Real> * stress;
akantu::Vector<akantu::Real> * strain;
akantu::Vector<akantu::Real> * damage;

#ifdef AKANTU_USE_IOHELPER
void paraviewInit(Dumper & dumper);
void paraviewDump(Dumper & dumper);
#endif

int main(int argc, char *argv[])
{
  akantu::UInt max_steps = 1000000;
  akantu::Real time_factor = 0.2;

  //  akantu::Real epot, ekin;

  akantu::Mesh mesh(spatial_dimension);
  akantu::MeshIOMSH mesh_io;
  mesh_io.read("bar2.msh", mesh);

  model = new akantu::SolidMechanicsModel(mesh);

  nb_nodes = model->getFEM().getMesh().getNbNodes();
  nb_element = model->getFEM().getMesh().getNbElement(type);


  /// model initialization
  model->initVectors();

  /// set vectors to 0
  memset(model->getForce().values,        0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getVelocity().values,     0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getAcceleration().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getDisplacement().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));


  model->initModel();
  model->readMaterials("material.dat");

  std::cout << model->getMaterial(0) << std::endl;
  std::cout << model->getMaterial(1) << std::endl;

  akantu::Vector<akantu::UInt> & elem_mat = model->getElementMaterial(type);
  for (akantu::UInt e = 0; e < nb_element; ++e) {
    akantu::Real barycenter[spatial_dimension];
    mesh.getBarycenter(e, type, barycenter, akantu::_not_ghost);
    if(barycenter[1] <= 0.75 && barycenter[1] >= 0.25)
      elem_mat(e, 0) = 0;
    else   elem_mat(e, 0) = 1;
  }

  model->initMaterials();
  model->assembleMassLumped();

  nb_quadrature_points = model->getFEM().getNbQuadraturePoints(type);

  stress = new akantu::Vector<akantu::Real>(nb_element * nb_quadrature_points,
					    spatial_dimension* spatial_dimension);
  strain = new akantu::Vector<akantu::Real>(nb_element * nb_quadrature_points,
					    spatial_dimension * spatial_dimension);
  damage = new akantu::Vector<akantu::Real>(nb_element * nb_quadrature_points, 1);
  memset(damage->values, 0, nb_quadrature_points * nb_element * sizeof(akantu::Real));

#ifdef AKANTU_USE_IOHELPER
  /// set to 0 only for the first paraview dump
  memset(model->getResidual().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  // memset(model->getMaterial(0).getStrain(type).values, 0,
  // 	 spatial_dimension*spatial_dimension*nb_element*sizeof(akantu::Real));
  // memset(model->getMaterial(0).getStress(type).values, 0,
  // 	 spatial_dimension*spatial_dimension*nb_element*sizeof(akantu::Real));
#endif //AKANTU_USE_IOHELPER


  /// boundary conditions
  akantu::Real eps = 1e-16;
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i] >= 9)
      model->getDisplacement().values[spatial_dimension*i] = (model->getFEM().getMesh().getNodes().values[spatial_dimension*i] - 9) / 10000. ;

    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i] <= eps)
	model->getBoundary().values[spatial_dimension*i] = true;

    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i + 1] <= eps ||
       model->getFEM().getMesh().getNodes().values[spatial_dimension*i + 1] >= 1 - eps ) {
      model->getBoundary().values[spatial_dimension*i + 1] = true;
    }
  }

  /// set the time step
  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);


#ifdef AKANTU_USE_IOHELPER
  /// initialize the paraview output
  DumperParaview dumper;
  paraviewInit(dumper);
#endif //AKANTU_USE_IOHELPER


#ifdef CHECK_STRESS
  std::ofstream outfile;
  outfile.open("stress");
#endif // CHECK_STRESS

  std::ofstream energy;
  energy.open("energy.csv");
  energy << "id,epot,ekin,tot" << std::endl;

  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    model->explicitPred();

    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    akantu::Real epot = model->getPotentialEnergy();
    akantu::Real ekin = model->getKineticEnergy();

    std::cerr << "passing step " << s << "/" << max_steps << std::endl;
    energy << s << "," << epot << "," << ekin << "," << epot + ekin
	   << std::endl;

#ifdef CHECK_STRESS
    /// search the position of the maximum of stress to determine the wave speed
    akantu::Real max_stress = std::numeric_limits<akantu::Real>::min();
    akantu::Real * stress = model->getMaterial(0).getStress(type).values;
    for (akantu::UInt i = 0; i < nb_element; ++i) {
      if(max_stress < stress[i*spatial_dimension*spatial_dimension]) {
	max_stress = stress[i*spatial_dimension*spatial_dimension];
      }
    }

    akantu::Real * coord    = model->getFEM().getMesh().getNodes().values;
    akantu::Real * disp_val = model->getDisplacement().values;
    akantu::UInt * conn     = model->getFEM().getMesh().getConnectivity(type).values;
    akantu::UInt nb_nodes_per_element = model->getFEM().getMesh().getNbNodesPerElement(type);
    akantu::Real * coords = new akantu::Real[spatial_dimension];
    akantu::Real min_x = std::numeric_limits<akantu::Real>::max();
    akantu::Real max_x = std::numeric_limits<akantu::Real>::min();

    akantu::Real stress_range = 5e7;
    for (akantu::UInt el = 0; el < nb_element; ++el) {
      if(stress[el*spatial_dimension*spatial_dimension] > max_stress - stress_range) {
	akantu::UInt el_offset  = el * nb_nodes_per_element;
	memset(coords, 0, spatial_dimension*sizeof(akantu::Real));
	for (akantu::UInt n = 0; n < nb_nodes_per_element; ++n) {
	  for (akantu::UInt i = 0; i < spatial_dimension; ++i) {
	    akantu::UInt node = conn[el_offset + n] * spatial_dimension;
	    coords[i] += (coord[node + i] + disp_val[node + i])
	      / ((akantu::Real) nb_nodes_per_element);
	  }
	}
	min_x = min_x < coords[0] ? min_x : coords[0];
	max_x = max_x > coords[0] ? max_x : coords[0];
      }
    }


    outfile << s << " " << .5 * (min_x + max_x) << " " << min_x << " " << max_x << " " << max_x - min_x << " " << max_stress << std::endl;

    delete [] coords;
#endif // CHECK_STRESS


#ifdef AKANTU_USE_IOHELPER
    if(s % 100 == 0) paraviewDump(dumper);
#endif //AKANTU_USE_IOHELPER
    if(s % 10 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;
  }

  energy.close();

#ifdef CHECK_STRESS
  outfile.close();
#endif // CHECK_STRESS

  delete model;

  akantu::finalize();

  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
/* Dumper vars                                                                */
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER
void paraviewInit(Dumper & dumper) {
  dumper.SetMode(TEXT);
  dumper.SetPoints(model->getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, "coordinates");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(type).values,
			 paraview_type, nb_element, C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().values,
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model->getVelocity().values,
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model->getResidual().values,
			  spatial_dimension, "force");
  dumper.AddNodeDataField(model->getMass().values,
			  1, "mass");
  dumper.AddNodeDataField(model->getForce().values,
			  spatial_dimension, "applied_force");

  akantu::Real * mat = new akantu::Real[nb_element * nb_quadrature_points];
  akantu::Vector<akantu::UInt> & elem_mat = model->getElementMaterial(type);
  for (akantu::UInt e = 0; e < nb_element; ++e) {
    for (akantu::UInt q = 0; q < nb_quadrature_points; ++q) {
      mat[e * nb_quadrature_points + q] = elem_mat(e, 0); 
    }
  }

  akantu::UInt offset = nb_quadrature_points * spatial_dimension * spatial_dimension;
  akantu::UInt nb_mat = model->getNbMaterials();
  for (akantu::UInt m = 0; m < nb_mat; ++m) {
    akantu::Material & material = model->getMaterial(m);
    const akantu::Vector<akantu::UInt> & elmat = material.getElementFilter(type);
    for (akantu::UInt e = 0; e < elmat.getSize(); ++e) {
      memcpy(stress->values + elmat(e, 0) * offset,
	     material.getStress(type).values + e * offset,
	     offset * sizeof(akantu::Real));
      memcpy(strain->values + elmat(e, 0) * offset,
	     material.getStrain(type).values + e * offset,
	     offset * sizeof(akantu::Real));
      if(m == 0)
	memcpy(damage->values + elmat(e, 0) * nb_quadrature_points,
	       dynamic_cast<akantu::MaterialDamage &>(material).getDamage(type).values + e * nb_quadrature_points,
	       nb_quadrature_points * sizeof(akantu::Real));

    }
  }

  dumper.AddElemDataField(mat, 1, "material");
  dumper.AddElemDataField(strain->values,
   			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(stress->values,
   			  spatial_dimension*spatial_dimension, "stress");
  dumper.AddElemDataField(damage->values,
   			  1, "damage");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("applied_force", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */
void paraviewDump(Dumper & dumper) {
  akantu::UInt offset = nb_quadrature_points * spatial_dimension * spatial_dimension;
  akantu::UInt nb_mat = model->getNbMaterials();
  for (akantu::UInt m = 0; m < nb_mat; ++m) {
    akantu::Material & material = model->getMaterial(m);
    const akantu::Vector<akantu::UInt> & elmat = material.getElementFilter(type);
    for (akantu::UInt e = 0; e < elmat.getSize(); ++e) {
      memcpy(stress->values + elmat(e, 0) * offset,
	     material.getStress(type).values + e * offset,
	     offset * sizeof(akantu::Real));
      memcpy(strain->values + elmat(e, 0) * offset,
	     material.getStrain(type).values + e * offset,
	     offset * sizeof(akantu::Real));
      if(m == 0)
	memcpy(damage->values + elmat(e, 0) * nb_quadrature_points,
	       dynamic_cast<akantu::MaterialDamage &>(material).getDamage(type).values + e * nb_quadrature_points,
	       nb_quadrature_points * sizeof(akantu::Real));

    }
  }

  dumper.Dump();
}
#endif
