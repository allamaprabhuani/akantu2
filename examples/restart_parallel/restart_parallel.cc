/**
 * @file   restart_parallel.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Sep  7 14:58:01 2011
 *
 * @brief  Example of using IOHelper by "hand"
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
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "static_communicator.hh"
#include "mesh_partition_scotch.hh"
#include "dof_synchronizer.hh"

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"

#endif //AKANTU_USE_IOHELPER

/* -------------------------------------------------------------------------- */

using namespace akantu;

#ifdef AKANTU_USE_IOHELPER
static void paraviewInit(Dumper & dumper,
			 const SolidMechanicsModel & model,
			 const ElementType & type,
			 const std::string & filename);
static void paraviewDump(Dumper & dumper);
static void checkpointInit(Dumper & dumper,
		    const SolidMechanicsModel & model,
		    const ElementType & type,
		    const std::string & filename);
static void checkpoint(Dumper & dumper,
		       const SolidMechanicsModel & model);
static void restart(const SolidMechanicsModel & model,
	     const ElementType & type,
	     const std::string & filename);
#endif

const UInt spatial_dimension = 3;

int main(int argc, char *argv[]) {
  debug::setDebugLevel(akantu::dblWarning);
  akantu::initialize(&argc, &argv);

  ElementType type = _tetrahedron_4;
  UInt max_steps = 100;

  Mesh mesh(spatial_dimension);

  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm->getNbProc();
  Int prank = comm->whoAmI();

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    MeshIOMSH mesh_io;
    mesh_io.read("cube.msh", mesh);
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }

  akantu::SolidMechanicsModel model(mesh);
  model.initParallel(partition);

  model.initModel();
  model.initVectors();
  model.initExplicit();
  model.readMaterials("material.dat");
  model.initMaterials();


  /// set vectors to 0
  model.getForce().clear();
  model.getVelocity().clear();
  model.getAcceleration().clear();
  model.getDisplacement().clear();

  /// boundary conditions
  Real eps = 1e-16;
  const akantu::Vector<akantu::Real> & pos = mesh.getNodes();
  akantu::Vector<akantu::Real> & force = model.getForce();
  akantu::Vector<bool> & boun = model.getBoundary();

  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    if(fabs(pos(i, 2) - 1) <= eps) force(i, 2) = -250;
    if(pos(i, 2) <= eps) boun(i, 2) = true;
    if(pos(i, 0) <= eps) boun(i, 0) = true;
  }



  model.assembleMassLumped();

  Real time_step = model.getStableTimeStep();
  if(prank == 0)
    std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step * .8);

  model.updateResidual();
#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  paraviewInit(dumper, model, type, "restart_parallel");
#endif

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    model.explicitPred();
    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

#ifdef AKANTU_USE_IOHELPER
    paraviewDump(dumper);

    if(s == (max_steps / 2)) {

      DumperParaview pv_dumper_restart;
      paraviewInit(pv_dumper_restart, model, type, "around_checkpoint_restart");

      DumperRestart dumper_restart;
      checkpointInit(dumper_restart, model, type, "restart");

      restart(model, type, "restart");
      model.updateResidual();
      paraviewDump(pv_dumper_restart);
    }
#endif //AKANTU_USE_IOHELPER
  }

  akantu::finalize();

  return EXIT_SUCCESS;
}



/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
template <ElementType type>
static UInt getIOHelperType() { AKANTU_DEBUG_TO_IMPLEMENT(); return 0; };

template <> UInt getIOHelperType<_segment_2>()      { return LINE1; }
template <> UInt getIOHelperType<_segment_3>()      { return LINE2; }
template <> UInt getIOHelperType<_triangle_3>()     { return TRIANGLE1; }
template <> UInt getIOHelperType<_triangle_6>()     { return TRIANGLE2; }
template <> UInt getIOHelperType<_quadrangle_4>()   { return QUAD1; }
template <> UInt getIOHelperType<_quadrangle_8>()   { return QUAD2; }
template <> UInt getIOHelperType<_tetrahedron_4>()  { return TETRA1; }
template <> UInt getIOHelperType<_tetrahedron_10>() { return TETRA2; }
template <> UInt getIOHelperType<_hexahedron_8>()   { return HEX1; }

static UInt getIOHelperType(ElementType type) {
  UInt ioh_type = 0;
#define GET_IOHELPER_TYPE(type)			\
  ioh_type = getIOHelperType<type>();

  AKANTU_BOOST_ELEMENT_SWITCH(GET_IOHELPER_TYPE);
#undef GET_IOHELPER_TYPE
  return ioh_type;
}

/* -------------------------------------------------------------------------- */
void paraviewInit(Dumper & dumper,
		  const SolidMechanicsModel & model,
		  const ElementType & type,
		  const std::string & filename) {
  const Mesh & mesh = model.getFEM().getMesh();
  UInt spatial_dimension = mesh.getSpatialDimension(type);
  UInt nb_nodes   = mesh.getNbNodes();
  UInt nb_element = mesh.getNbElement(type);

  std::stringstream filename_sstr; filename_sstr << filename << "_" << type;

  UInt whoami = StaticCommunicator::getStaticCommunicator()->whoAmI();
  UInt nproc  = StaticCommunicator::getStaticCommunicator()->getNbProc();

  dumper.SetMode(TEXT);
  dumper.SetParallelContext(whoami, nproc);

  dumper.SetPoints(mesh.getNodes().values,
		   spatial_dimension, nb_nodes, filename_sstr.str().c_str());
  dumper.SetConnectivity((int *)mesh.getConnectivity(type).values,
			 getIOHelperType(type), nb_element, C_MODE);

  dumper.AddNodeDataField(model.getDisplacement().values,
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model.getVelocity().values,
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model.getAcceleration().values,
			  spatial_dimension, "acceleration");
  dumper.AddNodeDataField(model.getMass().values,
			  spatial_dimension, "mass");
  dumper.AddNodeDataField(model.getResidual().values,
			  spatial_dimension, "force");
  dumper.AddNodeDataField(model.getForce().values,
			  spatial_dimension, "applied_force");

  dumper.AddElemDataField(model.getMaterial(0).getStrain(type).values,
   			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(model.getMaterial(0).getStrain(type).values,
   			  spatial_dimension*spatial_dimension, "stress");

  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("applied_force", 1);

  dumper.SetPrefix("paraview/");

  dumper.Init();
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */
void paraviewDump(Dumper & dumper) {
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

Vector<Real> checkpoint_displacements(0, spatial_dimension);
Vector<Real> checkpoint_velocity(0, spatial_dimension);
Vector<Real> checkpoint_acceleration(0, spatial_dimension);
Vector<Real> checkpoint_force(0, spatial_dimension);

/* -------------------------------------------------------------------------- */
void checkpointInit(Dumper & dumper,
		    const SolidMechanicsModel & model,
		    const ElementType & type,
		    const std::string & filename) {
  UInt whoami = StaticCommunicator::getStaticCommunicator()->whoAmI();
  UInt nproc  = StaticCommunicator::getStaticCommunicator()->getNbProc();

  Vector<Real> & displacements = model.getDisplacement();
  Vector<Real> & velocity      = model.getVelocity();
  Vector<Real> & acceleration  = model.getAcceleration();
  Vector<Real> & force         = model.getForce();

  DOFSynchronizer & dof_synchronizer = const_cast<DOFSynchronizer &>(model.getDOFSynchronizer());
  dof_synchronizer.initScatterGatherCommunicationScheme();

  if(whoami == 0){
    const Mesh & mesh = model.getFEM().getMesh();
    UInt nb_nodes   = mesh.getNbGlobalNodes();

    checkpoint_displacements.resize(nb_nodes);
    checkpoint_velocity	    .resize(nb_nodes);
    checkpoint_acceleration .resize(nb_nodes);
    checkpoint_force        .resize(nb_nodes);

    dof_synchronizer.gather(displacements, 0, &checkpoint_displacements);
    dof_synchronizer.gather(velocity     , 0, &checkpoint_velocity     );
    dof_synchronizer.gather(acceleration , 0, &checkpoint_acceleration );
    dof_synchronizer.gather(force        , 0, &checkpoint_force        );

    UInt nb_element = mesh.getNbElement(type);
    UInt spatial_dimension = mesh.getSpatialDimension(type);

    std::stringstream filename_sstr; filename_sstr << filename << "_" << type;

    dumper.SetMode(COMPRESSED);
    dumper.SetParallelContext(whoami, nproc);

    dumper.SetPoints(mesh.getNodes().values,
		     spatial_dimension, nb_nodes, filename_sstr.str().c_str());
    dumper.SetConnectivity((int *)mesh.getConnectivity(type).values,
			   getIOHelperType(type), nb_element, C_MODE);

    dumper.AddNodeDataField(checkpoint_displacements.storage(), spatial_dimension, "displacements");
    dumper.AddNodeDataField(checkpoint_velocity     .storage(), spatial_dimension, "velocity");
    dumper.AddNodeDataField(checkpoint_acceleration .storage(), spatial_dimension, "acceleration");
    dumper.AddNodeDataField(checkpoint_force        .storage(), spatial_dimension, "applied_force");
    dumper.SetPrefix("restart/");

    dumper.Init();
    dumper.Dump();
  } else {
    dof_synchronizer.gather(displacements, 0);
    dof_synchronizer.gather(velocity     , 0);
    dof_synchronizer.gather(acceleration , 0);
    dof_synchronizer.gather(force        , 0);
  }
}

/* -------------------------------------------------------------------------- */
void checkpoint(Dumper & dumper,
		const SolidMechanicsModel & model) {
  UInt whoami = StaticCommunicator::getStaticCommunicator()->whoAmI();

  DOFSynchronizer & dof_synchronizer = const_cast<DOFSynchronizer &>(model.getDOFSynchronizer());
  Vector<Real> & displacements = model.getDisplacement();
  Vector<Real> & velocity      = model.getVelocity();
  Vector<Real> & acceleration  = model.getAcceleration();
  Vector<Real> & force         = model.getForce();

  if(whoami == 0){
    dof_synchronizer.gather(displacements, 0, &checkpoint_displacements);
    dof_synchronizer.gather(velocity     , 0, &checkpoint_velocity     );
    dof_synchronizer.gather(acceleration , 0, &checkpoint_acceleration );
    dof_synchronizer.gather(force        , 0, &checkpoint_force        );

    dumper.Dump();
  } else {
    dof_synchronizer.gather(displacements, 0);
    dof_synchronizer.gather(velocity     , 0);
    dof_synchronizer.gather(acceleration , 0);
    dof_synchronizer.gather(force        , 0);
  }
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

void restart(const SolidMechanicsModel & model,
	     const ElementType & type,
	     const std::string & filename) {
  UInt whoami = StaticCommunicator::getStaticCommunicator()->whoAmI();
  UInt nproc  = StaticCommunicator::getStaticCommunicator()->getNbProc();

  DOFSynchronizer & dof_synchronizer = const_cast<DOFSynchronizer &>(model.getDOFSynchronizer());
  dof_synchronizer.initScatterGatherCommunicationScheme();

  Vector<Real> & displacements = model.getDisplacement();
  Vector<Real> & velocity      = model.getVelocity();
  Vector<Real> & acceleration  = model.getAcceleration();
  Vector<Real> & force         = model.getForce();

  if(whoami == 0){
    const Mesh & mesh = model.getFEM().getMesh();
    UInt nb_nodes   = mesh.getNbGlobalNodes();
    UInt spatial_dimension = mesh.getSpatialDimension(type);

    std::stringstream filename_sstr; filename_sstr << filename << "_" << type;

    ReaderRestart reader;

    reader.SetMode(COMPRESSED);
    reader.SetParallelContext(whoami, nproc);

    reader.SetPoints(filename_sstr.str().c_str());
    reader.SetConnectivity(getIOHelperType(type));

    reader.AddNodeDataField("displacements");
    reader.AddNodeDataField("velocity");
    reader.AddNodeDataField("acceleration");
    reader.AddNodeDataField("applied_force");
    reader.SetPrefix("restart/");

    reader.Init();
    reader.Read();

    Vector<Real> restart_tmp_vect(nb_nodes, spatial_dimension);

    displacements.clear();
    velocity.clear();
    acceleration.clear();
    force.clear();

    Real * tmp = reader.GetNodeDataField("displacements");
    std::copy(tmp, tmp + nb_nodes * spatial_dimension, restart_tmp_vect.storage());
    dof_synchronizer.scatter(displacements, 0, &restart_tmp_vect);

    tmp = reader.GetNodeDataField("velocity");
    std::copy(tmp, tmp + nb_nodes * spatial_dimension, restart_tmp_vect.storage());
    dof_synchronizer.scatter(velocity     , 0, &restart_tmp_vect);

    tmp = reader.GetNodeDataField("acceleration");
    std::copy(tmp, tmp + nb_nodes * spatial_dimension, restart_tmp_vect.storage());
    dof_synchronizer.scatter(acceleration , 0, &restart_tmp_vect);

    tmp = reader.GetNodeDataField("applied_force");
    std::copy(tmp, tmp + nb_nodes * spatial_dimension, restart_tmp_vect.storage());
    dof_synchronizer.scatter(force        , 0, &restart_tmp_vect);
  } else {
    displacements.clear();
    velocity.clear();
    acceleration.clear();
    force.clear();
    dof_synchronizer.scatter(displacements, 0);
    dof_synchronizer.scatter(velocity     , 0);
    dof_synchronizer.scatter(acceleration , 0);
    dof_synchronizer.scatter(force        , 0);
  }
}

/* -------------------------------------------------------------------------- */


#endif
