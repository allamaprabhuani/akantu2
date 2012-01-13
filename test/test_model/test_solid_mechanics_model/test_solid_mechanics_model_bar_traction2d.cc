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
#  include "io_helper.hh"


iohelper::ElemType paraview_type = iohelper::TRIANGLE2;
#endif //AKANTU_USE_IOHELPER

//#define CHECK_STRESS
akantu::ElementType type = akantu::_triangle_6;

akantu::SolidMechanicsModel * model;
akantu::UInt spatial_dimension = 2;
akantu::UInt nb_nodes;
akantu::UInt nb_element;
akantu::UInt nb_quadrature_points;

akantu::Vector<akantu::Real> * stress;
akantu::Vector<akantu::Real> * strain;
akantu::Vector<akantu::Real> * damage;

#ifdef AKANTU_USE_IOHELPER
static void paraviewInit(iohelper::Dumper & dumper);
static void paraviewDump(iohelper::Dumper & dumper);
#endif

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  akantu::UInt max_steps = 5000;
  akantu::Real time_factor = 0.8;

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
  model->getForce().clear();
  model->getVelocity().clear();
  model->getAcceleration().clear();
  model->getDisplacement().clear();

  model->initExplicit();
  model->initModel();
  model->readMaterials("material.dat");

  std::cout << model->getMaterial(0) << std::endl;

  model->initMaterials();
  model->assembleMassLumped();

  nb_quadrature_points = model->getFEM().getNbQuadraturePoints(type);

  stress = new akantu::Vector<akantu::Real>(nb_element * nb_quadrature_points,
					    spatial_dimension* spatial_dimension);
  strain = new akantu::Vector<akantu::Real>(nb_element * nb_quadrature_points,
					    spatial_dimension * spatial_dimension);

  /// boundary conditions
  akantu::Real eps = 1e-16;
  const akantu::Vector<akantu::Real> & pos = mesh.getNodes();
  akantu::Vector<akantu::Real> & disp = model->getDisplacement();
  akantu::Vector<bool> & boun = model->getBoundary();

  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if(pos(i, 0) >= 9.) disp(i, 0) = (pos(i, 0) - 9) / 100.;
    if(pos(i) <= eps)   boun(i, 0) = true;
    if(pos(i, 1) <= eps || pos(i, 1) >= 1 - eps ) boun(i, 1) = true;
  }

  /// set the time step
  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);


#ifdef AKANTU_USE_IOHELPER
  /// initialize the paraview output
  model->updateResidual();
  iohelper::DumperParaview dumper;
  paraviewInit(dumper);
#endif //AKANTU_USE_IOHELPER


#ifdef CHECK_STRESS
  std::ofstream outfile;
  outfile.open("stress");
#endif // CHECK_STRESS

  std::ofstream energy;
  energy.open("energy_bar_2d.csv");
  energy << "id,rtime,epot,ekin,tot" << std::endl;

  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    model->explicitPred();

    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    akantu::Real epot = model->getPotentialEnergy();
    akantu::Real ekin = model->getKineticEnergy();

    energy << s << "," << (s-1)*time_step << "," << epot << "," << ekin << "," << epot + ekin
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
    if(s % 100 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;
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
/* iohelper::Dumper vars                                                                */
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER
void paraviewInit(iohelper::Dumper & dumper) {
  dumper.SetMode(iohelper::TEXT);
  dumper.SetPoints(model->getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, "bar2d");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(type).values,
			 paraview_type, nb_element, iohelper::C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().values,
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model->getVelocity().values,
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model->getAcceleration().values,
			  spatial_dimension, "acceleration");
  dumper.AddNodeDataField(model->getResidual().values,
			  spatial_dimension, "force");
  dumper.AddNodeDataField(model->getMass().values,
			  spatial_dimension, "mass");
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
    }
  }

  dumper.AddElemDataField(mat, 1, "material");
  dumper.AddElemDataField(strain->values,
   			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(stress->values,
   			  spatial_dimension*spatial_dimension, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("applied_force", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */
void paraviewDump(iohelper::Dumper & dumper) {
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
    }
  }

  dumper.Dump();
}
#endif
