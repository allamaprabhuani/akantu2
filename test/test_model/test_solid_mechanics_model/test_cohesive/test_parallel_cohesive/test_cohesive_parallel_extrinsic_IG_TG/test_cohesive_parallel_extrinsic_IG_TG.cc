/**
 * @file   test_cohesive_extrinsic_IG_TG.cc
 *
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Thu April 18 13:31:00 2013
 *
 * @brief  Test for considering different cohesive properties for intergranular (IG) and
 * transgranular (TG) fractures in extrinsic cohesive elements
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
  const UInt max_steps = 1000;

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

  MultiGrainMaterialSelector material_selector(model, "TG_cohesive", "IG_cohesive");
  model.setMaterialSelector(material_selector);
  model.initFull("material.dat", SolidMechanicsModelCohesiveOptions(_explicit_lumped_mass, true, false));

  Real time_step = model.getStableTimeStep()*0.05;
  model.setTimeStep(time_step);
  //  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();


  /* ------------------------------------------------------------------------ */
  /* Facet part                                                               */
  /* ------------------------------------------------------------------------ */

  //  std::cout << mesh << std::endl;

  const Mesh & mesh_facets = model.getMeshFacets();

  Array<Real> & position = mesh.getNodes();

  Real * bary_facet = new Real[spatial_dimension];
  // first, the tag which shows grain ID should be read for each element

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end();
       ++gt) {
    GhostType gt_facet = *gt;
    Mesh::type_iterator first = mesh_facets.firstType(spatial_dimension - 1, gt_facet);
    Mesh::type_iterator last  = mesh_facets.lastType(spatial_dimension - 1, gt_facet);

    for(;first != last; ++first) {
      ElementType type_facet = *first;

      Array<bool> & facet_check = model.getFacetsCheck(type_facet);
      if (gt_facet == _not_ghost)
	facet_check.clear();

      UInt nb_facet = mesh_facets.getNbElement(type_facet, gt_facet);
      ////////////////////////////////////////////

      for (UInt f = 0; f < nb_facet; ++f) {
	mesh_facets.getBarycenter(f, type_facet, bary_facet, gt_facet);
	if ((bary_facet[1] < 0.1 && bary_facet[1] > -0.1) ||(bary_facet[0] < 0.1 && bary_facet[0] > -0.1)) {
	  if (gt_facet == _not_ghost)
	    facet_check(f) = true;
	}
      }
    }
  }
  delete[] bary_facet;

  // std::cout << nb_IG << " " << nb_TG << std::endl;

  //  model.initFacetFilter();

   //  for ( UInt i = 0; i < nb_facet; ++i){

  // }
  //  std::cout << mesh << std::endl;

  /* ------------------------------------------------------------------------ */
  /* End of facet part                                                        */
  /* ------------------------------------------------------------------------ */

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

    model.explicitPred();
    // nb_coh_elem = mesh.getNbElement (FEM::getCohesiveElementType(type_facet));
    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    model.dump();
    dumper.dump();
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

  // edis.close();
  // erev.close();

  //  mesh.write("mesh_final.msh");

  Real Ed = model.getEnergy("dissipated");

  Real Edt = 40;

  if(prank == 0)
    std::cout << Ed << " " << Edt << std::endl;

  if (Ed < Edt * 0.999 || Ed > Edt * 1.001 || std::isnan(Ed)) {
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
