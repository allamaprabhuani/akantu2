/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "contact_mechanics_model.hh"
#include "coupler_solid_contact.hh"
#include "non_linear_solver.hh"
#include "dumper_text.hh"
#include "dumper_variable.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {

  UInt max_steps = 50000;
  UInt imposing_steps = 5000;
  UInt damping_interval = 10;
  Real max_weight = 0.004;
  Real max_displacement = 0.01;
  Real damping_ratio = 0.9;
  std::string mesh_file = "hertz.msh";
  std::string material_file = "material_explicit.dat";
    
  const UInt spatial_dimension = 2;
  initialize(material_file, argc, argv);

  Real time_factor = 0.1; 
  
  Mesh mesh(spatial_dimension);
  mesh.read(mesh_file);
  
  CouplerSolidContact coupler(mesh);

  auto & solid   = coupler.getSolidMechanicsModel();
  auto & contact = coupler.getContactMechanicsModel();
  
  auto && material_selector = std::make_shared<MeshDataMaterialSelector<std::string>>(
			      "physical_names",solid);
  solid.setMaterialSelector(material_selector);

  solid.initFull(  _analysis_method = _explicit_lumped_mass);
  contact.initFull(_analysis_method = _explicit_dynamic_contact);

  auto &&surface_selector = std::make_shared<PhysicalSurfaceSelector>(mesh);
  contact.getContactDetector().setSurfaceSelector(surface_selector);
  
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "bottom");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "bottom");

  //Vector<Real> weight = {0, -max_weight};
  //solid.applyBC(BC::Neumann::FromSameDim(weight), "top");

  coupler.initFull(_analysis_method = _explicit_dynamic_contact);

  Real time_step = solid.getStableTimeStep();
  time_step *= time_factor;

  std::cout << "Time Step = " << time_step << "s (" << time_step
            << "s)" << std::endl;
  coupler.setTimeStep(time_step);
   
  coupler.setBaseName("explicit-dynamic-implicit");
  coupler.addDumpFieldVector("displacement");
  coupler.addDumpFieldVector("normals");
  coupler.addDumpFieldVector("tangents");
  coupler.addDumpFieldVector("contact_force");
  coupler.addDumpFieldVector("external_force");
  coupler.addDumpFieldVector("internal_force");
  coupler.addDumpField("gaps");
  coupler.addDumpField("areas");
  coupler.addDumpField("blocked_dofs");
  coupler.addDumpField("strain");
  coupler.addDumpField("stress");
  coupler.addDumpField("Von Mises stress");

  auto & velocity = solid.getVelocity();

  for (UInt s : arange(max_steps)) {

    std::cerr << "Step " << s << std::endl;  

    Real increment = max_displacement / Real(imposing_steps / 2);
    solid.applyBC(BC::Dirichlet::IncrementValue(-increment, _y), "top");
    
    coupler.solveStep();       
    
    if (s % damping_interval == 0) {
      for (auto & v : make_view(velocity)) 
	v *= damping_ratio;
    }

    if (s % 100 == 0) {
      coupler.dump();
    }
    
  }

  
  finalize();
  return EXIT_SUCCESS;
}

