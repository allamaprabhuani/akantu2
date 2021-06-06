/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "contact_mechanics_model.hh"
#include "coupler_solid_contact.hh"
#include "non_linear_solver.hh"
#include "surface_selector.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {

  
  const UInt spatial_dimension = 2;
  initialize("material.dat", argc, argv);

  Real time_step;
  Real time_factor = 0.1;
  UInt max_steps = 20000;
  Real max_displacement = 5e-3; 
    
  Mesh mesh(spatial_dimension);
  mesh.read("hertz.msh");
  
  CouplerSolidContact coupler(mesh);

  auto & solid   = coupler.getSolidMechanicsModel();
  auto & contact = coupler.getContactMechanicsModel();

  auto && selector = std::make_shared<MeshDataMaterialSelector<std::string>>(
		     "physical_names",solid);
  solid.setMaterialSelector(selector);
  
  coupler.initFull(  _analysis_method = _explicit_lumped_mass);
 
  auto && surface_selector = std::make_shared<PhysicalSurfaceSelector>(mesh);
  contact.getContactDetector().setSurfaceSelector(surface_selector);
 
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "fixed");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "fixed");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "loading");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "symmetry");
  
  time_step = solid.getStableTimeStep();
  time_step *= time_factor;
  std::cout << "Time Step = " << time_step << "s (" << time_step
            << "s)" << std::endl;
  coupler.setTimeStep(time_step);

  coupler.setBaseName("contact-explicit-dynamic");
  coupler.addDumpFieldVector("displacement");
  coupler.addDumpFieldVector("velocity");
  coupler.addDumpFieldVector("normals");
  coupler.addDumpFieldVector("contact_force");
  coupler.addDumpFieldVector("external_force");
  coupler.addDumpFieldVector("internal_force");
  coupler.addDumpField("gaps");
  coupler.addDumpField("areas");
  coupler.addDumpField("blocked_dofs");
  coupler.addDumpField("grad_u");
  coupler.addDumpField("stress");

  auto & velocity = solid.getVelocity();
  auto & gaps = contact.getGaps();
 
  Real damping_ratio = 0.99;
  auto increment = max_displacement/max_steps;
  
  for (auto i : arange(max_steps)) {
    
    solid.applyBC(BC::Dirichlet::IncrementValue(-increment, _y), "loading"); 

    coupler.solveStep();

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

  finalize();
  return EXIT_SUCCESS;
}

