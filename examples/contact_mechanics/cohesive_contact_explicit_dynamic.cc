/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model_cohesive.hh"
#include "contact_mechanics_model.hh"
#include "coupler_solid_cohesive_contact.hh"
#include "surface_selector.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]){

  const UInt spatial_dimension = 2;
  initialize("material-cohesive.dat", argc, argv);

  Real time_step{0.};
  Real time_factor = 0.1;
  UInt max_steps = 25000;
  Real max_displacement = 1e-3; 

  Mesh mesh(spatial_dimension);
  mesh.read("cohesive-contact.msh");

  CouplerSolidCohesiveContact coupler(mesh);

  auto & solid   = coupler.getSolidMechanicsModelCohesive();
  auto & contact = coupler.getContactMechanicsModel(); 

  auto && material_selector =
    std::make_shared<MeshDataMaterialCohesiveSelector>(solid);
  material_selector->setFallback(solid.getMaterialSelector());
  solid.setMaterialSelector(material_selector);
  
  auto && surface_selector =
    std::make_shared<CohesiveSurfaceSelector>(mesh);
  contact.getContactDetector().setSurfaceSelector(surface_selector);


  coupler.initFull(_analysis_method = _explicit_lumped_mass,
		   _is_extrinsic = true);

  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "sides");
  

  time_step = solid.getStableTimeStep();
  time_step *= time_factor;
  std::cout << "Time Step = " << time_step << "s (" << time_step
            << "s)" << std::endl;
  coupler.setTimeStep(time_step);

  coupler.setBaseName("cohesive-contact-explicit-dynamic");
  coupler.addDumpFieldVector("displacement");
  coupler.addDumpFieldVector("velocity");
  coupler.addDumpFieldVector("normals");
  coupler.addDumpField("blocked_dofs");
  coupler.addDumpField("grad_u");
  coupler.addDumpField("stress");
  coupler.addDumpField("gaps");
  coupler.addDumpField("areas");

  
  auto & velocity = solid.getVelocity();
  auto & gaps = contact.getGaps();
 
  Real damping_ratio = 0.99;
  auto increment = max_displacement/max_steps;

  for (auto i : arange(max_steps)) {
    
    solid.applyBC(BC::Dirichlet::IncrementValue( increment, _y), "loading"); 
    solid.applyBC(BC::Dirichlet::IncrementValue(-increment, _y), "fixed"); 

    coupler.solveStep();
    
    solid.checkCohesiveStress();
    
    // damping velocities only along the contacting zone    
    for(auto && tuple : zip(gaps,
			    make_view(velocity, spatial_dimension))){
      auto & gap = std::get<0>(tuple);
      auto & vel = std::get<1>(tuple);
      if(gap > 0) {
	vel *= damping_ratio;
      }
    }

    // dumping energies 
    if (i % 1000 == 0) {
      
      Real epot = solid.getEnergy("potential");
      Real ekin = solid.getEnergy("kinetic");

      std::cerr << i << "," << i * increment << "," << epot << "," << ekin << ","
		<< epot + ekin << "," << std::endl;
    }

    if (i % 1000 == 0) {
      coupler.dump();
    }
  }

  
  for (auto i : arange(max_steps)) {
    
    solid.applyBC(BC::Dirichlet::IncrementValue(-increment, _y), "loading"); 
    solid.applyBC(BC::Dirichlet::IncrementValue( increment, _y), "fixed"); 

    coupler.solveStep();
    
    solid.checkCohesiveStress();
    
    // damping velocities only along the contacting zone    
    for(auto && tuple : zip(gaps,
			    make_view(velocity, spatial_dimension))){
      auto & gap = std::get<0>(tuple);
      auto & vel = std::get<1>(tuple);
      if(gap > 0) {
	vel *= damping_ratio;
      }
    }

    // dumping energies 
    if (i % 1000 == 0) {
      
      Real epot = solid.getEnergy("potential");
      Real ekin = solid.getEnergy("kinetic");

      std::cerr << i << "," << i * increment << "," << epot << "," << ekin << ","
		<< epot + ekin << "," << std::endl;
    }

    if (i % 1000 == 0) {
      coupler.dump();
    }
  }


  
}
