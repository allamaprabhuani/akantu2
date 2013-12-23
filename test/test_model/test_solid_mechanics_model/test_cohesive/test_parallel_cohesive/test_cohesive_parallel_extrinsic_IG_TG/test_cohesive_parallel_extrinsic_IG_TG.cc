/**
 * @file   test_cohesive_extrinsic_IG_TG.cc
 *
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Thu April 18 13:31:00 2013
 *
 * @brief Test for considering different cohesive properties for
 * intergranular (IG) and transgranular (TG) fractures in extrinsic
 * cohesive elements
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
#include "material_cohesive_linear.hh"

#if defined(AKANTU_USE_IOHELPER)
#  include "dumper_paraview.hh"
#endif
/* -------------------------------------------------------------------------- */

using namespace akantu;

class MultiGrainMaterialSelector : public DefaultMaterialCohesiveSelector {
public:
  MultiGrainMaterialSelector(const SolidMechanicsModelCohesive & model, const ID & transgranular_id, const ID & intergranular_id) :
    DefaultMaterialCohesiveSelector(model),
    transgranular_id(transgranular_id),
    intergranular_id(intergranular_id),
    model(model),
    mesh(model.getMesh()),
    mesh_facets(model.getMeshFacets()),
    spatial_dimension(model.getSpatialDimension()),
    nb_IG(0), nb_TG(0) {
  }

  UInt operator()(const Element & element) {
    if(mesh_facets.getSpatialDimension(element.type) == (spatial_dimension - 1)) {
      const std::vector<Element> & element_to_subelement = mesh_facets.getElementToSubelement(element.type, element.ghost_type)(element.element);

      const Element & el1 = element_to_subelement[0];
      const Element & el2 = element_to_subelement[1];

      UInt grain_id1 = mesh.getData<UInt>("tag_0", el1.type, el1.ghost_type)(el1.element);
      if(el2 != ElementNull) {
	UInt grain_id2 = mesh.getData<UInt>("tag_0", el2.type, el2.ghost_type)(el2.element);
	if (grain_id1 == grain_id2){
	  //transgranular = 0 indicator
	  nb_TG++;
	  return model.getMaterialIndex(transgranular_id);
	} else  {
	  //intergranular = 1 indicator
	  nb_IG++;
	  return model.getMaterialIndex(intergranular_id);
	}
      } else {
	//transgranular = 0 indicator
	nb_TG++;
	return model.getMaterialIndex(transgranular_id);
      }
    } else {
      return DefaultMaterialCohesiveSelector::operator()(element);
    }
  }

private:
  ID transgranular_id, intergranular_id;
  const SolidMechanicsModelCohesive & model;
  const Mesh & mesh;
  const Mesh & mesh_facets;
  UInt spatial_dimension;

  UInt nb_IG;
  UInt nb_TG;
};


int main(int argc, char *argv[]) {
  initialize(argc, argv);

  debug::setDebugLevel(dblWarning);

  const UInt spatial_dimension = 2;
  const UInt max_steps = 500;

  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    mesh.read("square.msh");
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initParallel(partition, NULL, true);

  delete partition;

  MultiGrainMaterialSelector material_selector(model, "TG_cohesive", "IG_cohesive");
  model.setMaterialSelector(material_selector);
  model.initFull("material.dat", SolidMechanicsModelCohesiveOptions(_explicit_lumped_mass, true, false));

  Real time_step = model.getStableTimeStep()*0.1;
  model.setTimeStep(time_step);
  //  std::cout << "Time step: " << time_step << std::endl;

  //  std::cout << mesh << std::endl;

  Array<Real> & position = mesh.getNodes();
  Array<Real> & velocity = model.getVelocity();
  Array<bool> & boundary = model.getBoundary();
  Array<Real> & displacement = model.getDisplacement();
  //  const Array<Real> & residual = model.getResidual();

  UInt nb_nodes = mesh.getNbNodes();

  /// boundary conditions
  for (UInt n = 0; n < nb_nodes; ++n) {
    if (position(n, 1) > 0.99|| position(n, 1) < -0.99)
      boundary(n, 1) = true;

    if (position(n, 0) > 0.99 || position(n, 0) < -0.99)
      boundary(n, 0) = true;
  }

  model.synchronizeBoundaries();
  model.updateResidual();

  model.setBaseName("extrinsic");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("residual"    );
  model.addDumpField("stress");
  model.addDumpField("strain");
  model.addDumpField("partitions");

  DumperParaview dumper("cohesive_elements");
  dumper.registerMesh(mesh, spatial_dimension, _not_ghost, _ek_cohesive);
  DumperIOHelper::Field * cohesive_displacement =
    new DumperIOHelper::NodalField<Real>(model.getDisplacement());
  cohesive_displacement->setPadding(3);
  dumper.registerField("displacement", cohesive_displacement);

  dumper.dump();

  /// initial conditions
  Real loading_rate = 0.1;
  // bar_height  = 2
  Real VI = loading_rate * 2* 0.5;
  for (UInt n = 0; n < nb_nodes; ++n) {
    velocity(n, 1) = loading_rate * position(n, 1);
    velocity(n, 0) = loading_rate * position(n, 0);
  }

  // std::ofstream edis("edis.txt");
  // std::ofstream erev("erev.txt");

  //  Array<Real> & residual = model.getResidual();
  model.dump();
  //  const Array<Real> & stress = model.getMaterial(0).getStress(type);
  Real dispy = 0;
  // UInt nb_coh_elem = 0;

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {
    dispy += VI * time_step;
    /// update displacement on extreme nodes
    for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
      if (position(n, 1) > 0.99){
	displacement(n, 1) = dispy;
	velocity(n,1) = VI;}
      if (position(n, 1) < -0.99){
	displacement(n, 1) = -dispy;
	velocity(n,1) = -VI;}
      if (position(n, 0) > 0.99){
	displacement(n, 0) = dispy;
	velocity(n,0) = VI;}
      if (position(n, 0) < -0.99){
	displacement(n, 0) = -dispy;
	velocity(n,0) = -VI;}
    }

    model.checkCohesiveStress();

    model.solveStep();

    if(s % 10 == 0) {
      if(prank == 0)
	std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }

    // Real Ed = model.getEnergy("dissipated");

    // edis << s << " "
    //	 << Ed << std::endl;

    // erev << s << " "
    //	 << Er << std::endl;

  }

  model.dump();
  dumper.dump();

  // edis.close();
  // erev.close();

  //  mesh.write("mesh_final.msh");

  Real Ed = model.getEnergy("dissipated");

  Real Edt = 40;

  if(prank == 0)
    std::cout << Ed << " " << Edt << std::endl;

  if (Ed < Edt * 0.99 || Ed > Edt * 1.01 || std::isnan(Ed)) {
    if(prank == 0)
      std::cout << "The dissipated energy is incorrect" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  // for (UInt n = 0; n < position.getSize(); ++n) {
  //   for (UInt s = 0; s < spatial_dimension; ++s) {
  //     position(n, s) += displacement(n, s);
  //   }
  // }


  finalize();

  if(prank == 0)
    std::cout << "OK: test_cohesive_extrinsic_IG_TG was passed!" << std::endl;
  return EXIT_SUCCESS;
}
