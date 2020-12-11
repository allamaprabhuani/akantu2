/* -------------------------------------------------------------------------- */
#include "contact_mechanics_model.hh"
#include "coupler_solid_contact.hh"
#include "non_linear_solver.hh"
#include "solid_mechanics_model.hh"
#include "surface_selector.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {

  UInt max_steps        = 20000;
  Real max_weight       = 0.02;
  UInt damping_interval = 10;
  Real damping_ratio    = 0.9;

  std::string mesh_file     = "flat_on_flat.msh";
  std::string material_file = "material.dat";

  const UInt spatial_dimension = 2;

  initialize(material_file, argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read(mesh_file);

  CouplerSolidContact coupler(mesh);

  auto & solid = coupler.getSolidMechanicsModel();
  auto & contact = coupler.getContactMechanicsModel();

  auto && material_selector =
      std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names",
                                                              solid);
  solid.setMaterialSelector(material_selector);

  solid.initFull(_analysis_method   = _explicit_lumped_mass);
  contact.initFull(_analysis_method = _explicit_dynamic_contact);

  auto && surface_selector = std::make_shared<PhysicalSurfaceSelector>(mesh);
  contact.getContactDetector().setSurfaceSelector(surface_selector);

  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "bottom");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "bottom");

  Vector<Real> weight = {0, -max_weight};
  solid.applyBC(BC::Neumann::FromSameDim(weight), "top");

  coupler.initFull(_analysis_method = _explicit_dynamic_contact);

  Real time_step = solid.getStableTimeStep();
  time_step *= 0.1;

  coupler.setTimeStep(time_step);

  coupler.setBaseName("flat-on-flat");
  coupler.addDumpFieldVector("displacement");
  coupler.addDumpFieldVector("normals");
  coupler.addDumpFieldVector("contact_force");
  coupler.addDumpFieldVector("external_force");
  coupler.addDumpFieldVector("internal_force");
  coupler.addDumpField("gaps");
  coupler.addDumpField("areas");
  coupler.addDumpField("blocked_dofs");
  coupler.addDumpField("strain");
  coupler.addDumpField("stress");

  auto & velocity = solid.getVelocity();

  for (UInt s : arange(max_steps)) {

    coupler.solveStep();

    if (s % damping_interval == 0) {
      for (auto & v : make_view(velocity))
        v *= damping_ratio;
    }

    if (s % 100 == 0) {
      coupler.dump();
    }
  }

  const ElementType element_type =  _quad_4;
  const Array<Real> & stress_vect = model.getMaterial("top_body").getStress(element_type);
  
  auto stress_it = stress_vect.begin(spatial_dimension, spatial_dimension);
  auto stress_end = stress_vect.end(spatial_dimension, spatial_dimension);

  Real stress_tolerance = 1e-13;

  Matrix<Real> presc_stress;
  
  for (; stress_it != stress_end; ++stress_it) {
    const auto & stress = *stress_it;
    Matrix<Real> diff(spatial_dimension, spatial_dimension);

    diff = stress;
    diff -= presc_stress;
    Real stress_error = diff.norm<L_inf>() / stress.norm<L_inf>();

    if (stress_error > stress_tolerance) {
      std::cerr << "stress error: " << stress_error << " > " << stress_tolerance
                << std::endl;
      std::cerr << "stress: " << stress << std::endl
                << "prescribed stress: " << presc_stress << std::endl;
      return EXIT_FAILURE;
    }
  }
  
  finalize();
  return EXIT_SUCCESS;
}
