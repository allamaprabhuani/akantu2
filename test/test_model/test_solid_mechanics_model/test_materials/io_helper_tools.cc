/**
 * @file   io_helper_tools.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Sep 30 11:13:01 2011
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
#include "io_helper_tools.hh"

#include "aka_common.hh"
#include "mesh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "static_communicator.hh"
// #include "dof_synchronizer.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

/* ------------------------------------------------------------------------ */
template <ElementType type>
static iohelper::ElemType getIOHelperType() { AKANTU_DEBUG_TO_IMPLEMENT(); return iohelper::MAX_ELEM_TYPE; };

template <> iohelper::ElemType getIOHelperType<_segment_2>()      { return iohelper::LINE1; }
template <> iohelper::ElemType getIOHelperType<_segment_3>()      { return iohelper::LINE2; }
template <> iohelper::ElemType getIOHelperType<_triangle_3>()     { return iohelper::TRIANGLE1; }
template <> iohelper::ElemType getIOHelperType<_triangle_6>()     { return iohelper::TRIANGLE2; }
template <> iohelper::ElemType getIOHelperType<_quadrangle_4>()   { return iohelper::QUAD1; }
template <> iohelper::ElemType getIOHelperType<_quadrangle_8>()   { return iohelper::QUAD2; }
template <> iohelper::ElemType getIOHelperType<_tetrahedron_4>()  { return iohelper::TETRA1; }
template <> iohelper::ElemType getIOHelperType<_tetrahedron_10>() { return iohelper::TETRA2; }
template <> iohelper::ElemType getIOHelperType<_hexahedron_8>()   { return iohelper::HEX1; }

iohelper::ElemType getIOHelperType(ElementType type) {
  iohelper::ElemType ioh_type = iohelper::MAX_ELEM_TYPE;
#define GET_IOHELPER_TYPE(type)			\
  ioh_type = getIOHelperType<type>();

  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_IOHELPER_TYPE);
#undef GET_IOHELPER_TYPE
  return ioh_type;
}

/* ------------------------------------------------------------------------ */
void paraviewInit(iohelper::Dumper & dumper,
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

  dumper.SetMode(iohelper::TEXT);
  dumper.SetParallelContext(whoami, nproc);

  dumper.SetPoints(mesh.getNodes().values,
		   spatial_dimension, nb_nodes, filename_sstr.str().c_str());
  dumper.SetConnectivity((int *)mesh.getConnectivity(type).values,
			 getIOHelperType(type), nb_element, iohelper::C_MODE);

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
  dumper.AddElemDataField(model.getMaterial(0).getStress(type).values,
			  spatial_dimension*spatial_dimension, "stress");


  if(dynamic_cast<const MaterialDamage *>(&model.getMaterial(0))) {
    const MaterialDamage & mat = dynamic_cast<const MaterialDamage &>(model.getMaterial(0));
    Real * dam = mat.getDamage(type).storage();
    dumper.AddElemDataField(dam, 1, "damage");
  }

  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("applied_force", 1);

  dumper.SetPrefix("paraview/");

  dumper.Init();
  dumper.Dump();
}

/* ------------------------------------------------------------------------ */
void paraviewDump(iohelper::Dumper & dumper) {
  dumper.Dump();
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

// Vector<Real> checkpoint_displacements(0, spatial_dimension);
// Vector<Real> checkpoint_velocity(0, spatial_dimension);
// Vector<Real> checkpoint_acceleration(0, spatial_dimension);
// Vector<Real> checkpoint_force(0, spatial_dimension);

// /* ------------------------------------------------------------------------ */
// void checkpointInit(iohelper::Dumper & dumper,
// 		    const SolidMechanicsModel & model,
// 		    const ElementType & type,
// 		    const std::string & filename) {
//   UInt whoami = StaticCommunicator::getStaticCommunicator()->whoAmI();
//   UInt nproc  = StaticCommunicator::getStaticCommunicator()->getNbProc();

//   Vector<Real> & displacements = model.getDisplacement();
//   Vector<Real> & velocity      = model.getVelocity();
//   Vector<Real> & acceleration  = model.getAcceleration();
//   Vector<Real> & force         = model.getForce();

//   DOFSynchronizer & dof_synchronizer = const_cast<DOFSynchronizer &>(model.getDOFSynchronizer());
//   dof_synchronizer.initScatterGatherCommunicationScheme();

//   if(whoami == 0){
//     const Mesh & mesh = model.getFEM().getMesh();
//     UInt nb_nodes   = mesh.getNbGlobalNodes();

//     checkpoint_displacements.resize(nb_nodes);
//     checkpoint_velocity	    .resize(nb_nodes);
//     checkpoint_acceleration .resize(nb_nodes);
//     checkpoint_force        .resize(nb_nodes);

//     dof_synchronizer.gather(displacements, 0, &checkpoint_displacements);
//     dof_synchronizer.gather(velocity     , 0, &checkpoint_velocity     );
//     dof_synchronizer.gather(acceleration , 0, &checkpoint_acceleration );
//     dof_synchronizer.gather(force        , 0, &checkpoint_force        );

//     UInt nb_element = mesh.getNbElement(type);
//     UInt spatial_dimension = mesh.getSpatialDimension(type);

//     std::stringstream filename_sstr; filename_sstr << filename << "_" << type;

//     dumper.SetMode(iohelper::COMPRESSED);
//     dumper.SetParallelContext(whoami, nproc);

//     dumper.SetPoints(mesh.getNodes().values,
// 		     spatial_dimension, nb_nodes, filename_sstr.str().c_str());
//     dumper.SetConnectivity((int *)mesh.getConnectivity(type).values,
// 			   getIOHelperType(type), nb_element, iohelper::C_MODE);

//     dumper.AddNodeDataField(checkpoint_displacements.storage(), spatial_dimension, "displacements");
//     dumper.AddNodeDataField(checkpoint_velocity     .storage(), spatial_dimension, "velocity");
//     dumper.AddNodeDataField(checkpoint_acceleration .storage(), spatial_dimension, "acceleration");
//     dumper.AddNodeDataField(checkpoint_force        .storage(), spatial_dimension, "applied_force");
//     dumper.SetPrefix("restart/");

//     dumper.Init();
//     dumper.Dump();
//   } else {
//     dof_synchronizer.gather(displacements, 0);
//     dof_synchronizer.gather(velocity     , 0);
//     dof_synchronizer.gather(acceleration , 0);
//     dof_synchronizer.gather(force        , 0);
//   }
// }

// /* ------------------------------------------------------------------------ */
// void checkpoint(iohelper::Dumper & dumper,
// 		const SolidMechanicsModel & model) {
//   UInt whoami = StaticCommunicator::getStaticCommunicator()->whoAmI();

//   DOFSynchronizer & dof_synchronizer = const_cast<DOFSynchronizer &>(model.getDOFSynchronizer());
//   Vector<Real> & displacements = model.getDisplacement();
//   Vector<Real> & velocity      = model.getVelocity();
//   Vector<Real> & acceleration  = model.getAcceleration();
//   Vector<Real> & force         = model.getForce();

//   if(whoami == 0){
//     dof_synchronizer.gather(displacements, 0, &checkpoint_displacements);
//     dof_synchronizer.gather(velocity     , 0, &checkpoint_velocity     );
//     dof_synchronizer.gather(acceleration , 0, &checkpoint_acceleration );
//     dof_synchronizer.gather(force        , 0, &checkpoint_force        );

//     dumper.Dump();
//   } else {
//     dof_synchronizer.gather(displacements, 0);
//     dof_synchronizer.gather(velocity     , 0);
//     dof_synchronizer.gather(acceleration , 0);
//     dof_synchronizer.gather(force        , 0);
//   }
// }

// /* ------------------------------------------------------------------------ */
// /* ------------------------------------------------------------------------ */

// void restart(const SolidMechanicsModel & model,
// 	     const ElementType & type,
// 	     const std::string & filename) {
//   UInt whoami = StaticCommunicator::getStaticCommunicator()->whoAmI();
//   UInt nproc  = StaticCommunicator::getStaticCommunicator()->getNbProc();

//   DOFSynchronizer & dof_synchronizer = const_cast<DOFSynchronizer &>(model.getDOFSynchronizer());
//   dof_synchronizer.initScatterGatherCommunicationScheme();

//   Vector<Real> & displacements = model.getDisplacement();
//   Vector<Real> & velocity      = model.getVelocity();
//   Vector<Real> & acceleration  = model.getAcceleration();
//   Vector<Real> & force         = model.getForce();

//   if(whoami == 0){
//     const Mesh & mesh = model.getFEM().getMesh();
//     UInt nb_nodes   = mesh.getNbGlobalNodes();
//     UInt spatial_dimension = mesh.getSpatialDimension(type);

//     std::stringstream filename_sstr; filename_sstr << filename << "_" << type;

//     iohelper::ReaderRestart reader;

//     reader.SetMode(iohelper::COMPRESSED);
//     reader.SetParallelContext(whoami, nproc);

//     reader.SetPoints(filename_sstr.str().c_str());
//     reader.SetConnectivity(getIOHelperType(type));

//     reader.AddNodeDataField("displacements");
//     reader.AddNodeDataField("velocity");
//     reader.AddNodeDataField("acceleration");
//     reader.AddNodeDataField("applied_force");
//     reader.SetPrefix("restart/");

//     reader.Init();
//     reader.Read();

//     Vector<Real> restart_tmp_vect(nb_nodes, spatial_dimension);

//     displacements.clear();
//     velocity.clear();
//     acceleration.clear();
//     force.clear();

//     Real * tmp = reader.GetNodeDataField("displacements");
//     std::copy(tmp, tmp + nb_nodes * spatial_dimension, restart_tmp_vect.storage());
//     dof_synchronizer.scatter(displacements, 0, &restart_tmp_vect);

//     tmp = reader.GetNodeDataField("velocity");
//     std::copy(tmp, tmp + nb_nodes * spatial_dimension, restart_tmp_vect.storage());
//     dof_synchronizer.scatter(velocity     , 0, &restart_tmp_vect);

//     tmp = reader.GetNodeDataField("acceleration");
//     std::copy(tmp, tmp + nb_nodes * spatial_dimension, restart_tmp_vect.storage());
//     dof_synchronizer.scatter(acceleration , 0, &restart_tmp_vect);

//     tmp = reader.GetNodeDataField("applied_force");
//     std::copy(tmp, tmp + nb_nodes * spatial_dimension, restart_tmp_vect.storage());
//     dof_synchronizer.scatter(force        , 0, &restart_tmp_vect);
//   } else {
//     displacements.clear();
//     velocity.clear();
//     acceleration.clear();
//     force.clear();
//     dof_synchronizer.scatter(displacements, 0);
//     dof_synchronizer.scatter(velocity     , 0);
//     dof_synchronizer.scatter(acceleration , 0);
//     dof_synchronizer.scatter(force        , 0);
//   }
// }
/* ------------------------------------------------------------------------ */
