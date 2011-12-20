/**
 * @file   test_material_damage_non_local.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Sep  8 17:31:42 2011
 *
 * @brief  
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "fem.hh"
#include "local_material_damage.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"
using namespace iohelper;
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

akantu::Real eps = 1e-10;

static void trac(__attribute__ ((unused)) Real * position,
		 double * stress,
		 __attribute__ ((unused)) Real * normal,
		 __attribute__ ((unused)) UInt surface_id){
  memset(stress, 0, sizeof(Real)*4);
  if (fabs(position[0] - 10) < eps){
    stress[0] = 3e6;
    stress[3] = 3e6;
  }
}

int main(int argc, char *argv[])
{
  debug::setDebugLevel(dblWarning);

  akantu::initialize(argc, argv);
  UInt max_steps = 40000;

  Real bar_height = 4.;

  const UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  MeshIOMSH mesh_io;

  mesh_io.read("barre_trou.msh", mesh);

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initVectors();
  UInt nb_nodes = model.getFEM().getMesh().getNbNodes();


  model.initExplicit();
  model.initModel();
  model.readMaterials("material_damage_non_local.dat");
  model.initMaterials();

  /// set vectors to 0
  model.getForce().clear();
  model.getVelocity().clear();
  model.getAcceleration().clear();
  model.getDisplacement().clear();



  Real time_step = model.getStableTimeStep();
  model.setTimeStep(time_step/10.);

  model.assembleMassLumped();

  std::cout << model << std::endl;

  /// Dirichlet boundary conditions
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {

    if(model.getFEM().getMesh().getNodes().values[spatial_dimension*i] <= eps)
	model.getBoundary().values[spatial_dimension*i] = true;

    if(model.getFEM().getMesh().getNodes().values[spatial_dimension*i + 1] <= eps ||
       model.getFEM().getMesh().getNodes().values[spatial_dimension*i + 1] >= bar_height - eps ) {
      model.getBoundary().values[spatial_dimension*i + 1] = true;
    }
  }


  FEM & fem_boundary = model.getFEMBoundary();
  fem_boundary.initShapeFunctions();
  fem_boundary.computeNormalsOnControlPoints();
  model.computeForcesFromFunction(trac, akantu::_bft_stress);

  MaterialDamage & mat = dynamic_cast<MaterialDamage &>((model.getMaterial(0)));

#ifdef AKANTU_USE_IOHELPER
  model.updateResidual();

  DumperParaview dumper;
  dumper.SetMode(BASE64);

  dumper.SetPoints(model.getFEM().getMesh().getNodes().values, 2, nb_nodes, "coordinates_damage_nl");
  dumper.SetConnectivity((int *)model.getFEM().getMesh().getConnectivity(_triangle_6).values,
			 TRIANGLE2, model.getFEM().getMesh().getNbElement(_triangle_6), C_MODE);
  dumper.AddNodeDataField(model.getDisplacement().values, 2, "displacements");
  dumper.AddNodeDataField(model.getVelocity().values, 2, "velocity");
  dumper.AddNodeDataField(model.getForce().values, 2, "force");
  dumper.AddNodeDataField(model.getMass().values, 1, "Mass");
  dumper.AddNodeDataField(model.getResidual().values, 2, "residual");
  dumper.AddElemDataField(mat.getStrain(_triangle_6).values, 4, "strain");
  dumper.AddElemDataField(mat.getStress(_triangle_6).values, 4, "stress");

  Real * dam = mat.getDamage(_triangle_6).values;
  dumper.AddElemDataField(dam, 1, "damage");

  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("force", 1);
  dumper.SetEmbeddedValue("residual", 1);
  dumper.SetEmbeddedValue("velocity", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  //  mat.savePairs("cl_pairs");

  for(UInt s = 0; s < max_steps; ++s) {
    model.explicitPred();

    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    if(s % 100 == 0) std::cout << "Step " << s+1 << "/" << max_steps <<std::endl;

#ifdef AKANTU_USE_IOHELPER
    if(s % 100 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  }

  akantu::finalize();
  return EXIT_SUCCESS;
}
