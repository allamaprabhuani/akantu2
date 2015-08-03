/**
 * @file   test_local_material.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Clement Roux-Langlois <clement.roux@epfl.ch>
 *
 * @date creation: Fri Nov 26 2010
 * @date last modification: Fri Sep 19 2014
 *
 * @brief  test of the class SolidMechanicsModel
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "local_material_damage.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[])
{
  akantu::initialize("material.dat", argc, argv);
  UInt max_steps = 1100;
  Real epot, ekin;

  const UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  mesh.read("mesh_section_gap.msh");
  mesh.createGroupsFromMeshData<std::string>("physical_names");

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull(SolidMechanicsModelOptions(_explicit_lumped_mass, true));
  model.registerNewCustomMaterials<LocalMaterialDamage>("local_damage");
  model.initMaterials();

  Real time_step = model.getStableTimeStep();
  model.setTimeStep(time_step/2.5);

  model.assembleMassLumped();

  std::cout << model << std::endl;

  /// Dirichlet boundary conditions
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "Fixed");
  // model.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "Fixed");

  // Boundary condition (Neumann)
  Matrix<Real> stress(2,2);
  stress.eye(7e5);
  model.applyBC(BC::Neumann::FromHigherDim(stress), "Traction");

  /*model.setBaseName("damage_local");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("force"       );
  model.addDumpField("residual"    );
  model.addDumpField("damage"      );
  model.addDumpField("stress"      );
  model.addDumpField("strain"      );
  model.dump();*/

  for(UInt s = 0; s < max_steps; ++s) {
    if(s < 100){
    // Boundary condition (Neumann)
      stress.eye(7e5);
    model.applyBC(BC::Neumann::FromHigherDim(stress), "Traction");
    }

    model.solveStep();

    /*epot = model.getPotentialEnergy();
    ekin = model.getKineticEnergy();

    if(s % 10 == 0) std::cout << s << " " << epot << " " << ekin << " " << epot + ekin
	      << std::endl;

    if(s % 10 == 0) std::cout << "Step " << s+1 << "/" << max_steps <<std::endl;
    if(s % 10 == 0) model.dump();*/
  }

  const Vector<Real> & lower_bounds = mesh.getLowerBounds();
  const Vector<Real> & upper_bounds = mesh.getUpperBounds();
  Real L = upper_bounds(0) - lower_bounds(0);

  const ElementTypeMapArray<UInt> & filter =  model.getMaterial(0).getElementFilter();
  ElementTypeMapArray<UInt>::type_iterator it = filter.firstType(spatial_dimension);
  ElementTypeMapArray<UInt>::type_iterator end = filter.lastType(spatial_dimension);
  Vector<Real> barycenter(spatial_dimension);
  bool is_fully_damaged = false;
  for(; it != end; ++it) {
    UInt nb_elem = mesh.getNbElement(*it);
    const UInt nb_gp = model.getFEEngine().getNbQuadraturePoints(*it);
    Array<Real> & material_damage_array = model.getMaterial(0).getArray("damage", *it);
    UInt cpt = 0;
    for(UInt nel = 0; nel < nb_elem ; ++nel){
      if (material_damage_array(cpt,0) > 0.9){
	is_fully_damaged = true;
	mesh.getBarycenter(nel,*it,barycenter.storage());
	if( (std::abs(barycenter(0)-(L/2)) < (L/10) ) ) {
	  return EXIT_FAILURE;
	}
      }
      cpt += nb_gp;
    }
  }
  if(!is_fully_damaged)
    return EXIT_FAILURE;

  akantu::finalize();
  return EXIT_SUCCESS;
}
