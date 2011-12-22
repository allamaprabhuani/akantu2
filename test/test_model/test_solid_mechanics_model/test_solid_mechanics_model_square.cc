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
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "fem.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"
using iohelper::ElemType;
using iohelper::DumperParaview;
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

static void trac(__attribute__ ((unused)) double * position,
		 double * traction,
		 __attribute__ ((unused)) Real * normal,
		 __attribute__ ((unused)) UInt surface_id){
  memset(traction,0,sizeof(Real)*4);
  traction[0] = 1000;
  traction[3] = 1000;

  // if(fabs(position[0])< 1e-4){
  //   traction[0] = -position[1];
  // }

}

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  UInt max_steps = 1000;
  Real epot, ekin;

#ifdef AKANTU_USE_IOHELPER
  ElementType type  = _triangle_3;
  iohelper::ElemType para_type = iohelper::TRIANGLE1;
#endif //AKANTU_USE_IOHELPER

  Mesh mesh(2);
  MeshIOMSH mesh_io;
  mesh_io.read("square.msh", mesh);

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initVectors();
  UInt nb_nodes = model.getFEM().getMesh().getNbNodes();
  memset(model.getForce().values,        0, 2*nb_nodes*sizeof(Real));
  memset(model.getVelocity().values,     0, 2*nb_nodes*sizeof(Real));
  memset(model.getAcceleration().values, 0, 2*nb_nodes*sizeof(Real));
  memset(model.getDisplacement().values, 0, 2*nb_nodes*sizeof(Real));
  memset(model.getResidual().values,     0, 2*nb_nodes*sizeof(Real));
  memset(model.getMass().values,     1, nb_nodes*sizeof(Real));

  model.initExplicit();

  model.initModel();
  model.readMaterials("material.dat");
  model.initMaterials();


  Real time_step = model.getStableTimeStep();
  model.setTimeStep(time_step/10.);

  model.assembleMassLumped();

  std::cout << model << std::endl;

  /// boundary conditions
  // Real eps = 1e-16;
  // for (UInt i = 0; i < nb_nodes; ++i) {
  //   model.getDisplacement().values[2*i] = model.getFEM().getMesh().getNodes().values[2*i] / 100.;

  //   if(model.getFEM().getMesh().getNodes().values[2*i] <= eps) {
  //     model.getBoundary().values[2*i    ] = true;
  //     if(model.getFEM().getMesh().getNodes().values[2*i + 1] <= eps)
  // 	model.getBoundary().values[2*i + 1] = true;
  //   }
  //   if(model.getFEM().getMesh().getNodes().values[2*i + 1] <= eps) {
  //     model.getBoundary().values[2*i + 1] = true;
  //   }

  // }

  FEM & fem_boundary = model.getFEMBoundary();
  fem_boundary.initShapeFunctions();
  fem_boundary.computeNormalsOnControlPoints();
  model.computeForcesFromFunction(trac, akantu::_bft_stress);


#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperParaview dumper;
  dumper.SetMode(iohelper::BASE64);

  dumper.SetPoints(model.getFEM().getMesh().getNodes().values, 2, nb_nodes, "coordinates");
  dumper.SetConnectivity((int *)model.getFEM().getMesh().getConnectivity(type).values,
			 para_type, model.getFEM().getMesh().getNbElement(type), iohelper::C_MODE);
  dumper.AddNodeDataField(model.getDisplacement().values, 2, "displacements");
  dumper.AddNodeDataField(model.getVelocity().values, 2, "velocity");
  dumper.AddNodeDataField(model.getForce().values, 2, "force");
  dumper.AddNodeDataField(model.getMass().values, 1, "Mass");
  dumper.AddNodeDataField(model.getResidual().values, 2, "residual");
  dumper.AddElemDataField(model.getMaterial(0).getStrain(type).values, 4, "strain");
  dumper.AddElemDataField(model.getMaterial(0).getStress(type).values, 4, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("force", 1);
  dumper.SetEmbeddedValue("residual", 1);
  dumper.SetEmbeddedValue("velocity", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  std::ofstream energy;
  energy.open("energy.csv");
  energy << "id,epot,ekin,tot" << std::endl;

  for(UInt s = 0; s < max_steps; ++s) {
    model.explicitPred();

    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    epot = model.getPotentialEnergy();
    ekin = model.getKineticEnergy();

    std::cerr << "passing step " << s << "/" << max_steps << std::endl;
    energy << s << "," << epot << "," << ekin << "," << epot + ekin
	   << std::endl;

#ifdef AKANTU_USE_IOHELPER
    if(s % 100 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  }

  energy.close();

  finalize();

  return EXIT_SUCCESS;
}



