option(AKANTU_CONTACT "Use Contact package of Akantu" OFF)
set(CONTACT_FILES
#cc files

  model/solid_mechanics/contact.cc
  model/solid_mechanics/contact_search.cc
  model/solid_mechanics/contact_neighbor_structure.cc
  model/solid_mechanics/contact/contact_2d_explicit.cc
  model/solid_mechanics/contact/contact_search_2d_explicit.cc
  model/solid_mechanics/contact/regular_grid_neighbor_structure.cc
  model/solid_mechanics/contact/contact_search_explicit.cc
  model/solid_mechanics/contact/contact_3d_explicit.cc
  model/solid_mechanics/contact/grid_2d_neighbor_structure.cc
  model/solid_mechanics/contact/contact_rigid.cc
  model/solid_mechanics/contact/friction_coefficient.cc
  model/solid_mechanics/contact/friction_coefficient/unique_constant_fric_coef.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_weakening_coulomb.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/velocity_weakening_exponential.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/rice_kuwano.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/rice_kuwano_modified.cc
  model/solid_mechanics/contact/friction_coefficient/simplified_dieterich_fric_coef.cc
  model/solid_mechanics/contact/friction_coefficient/simplified_dieterich_fric_coef/ruina_slowness_fric_coef.cc
  model/solid_mechanics/contact/friction_coefficient/historic_velocity_fric_coef.cc

# include files

  model/solid_mechanics/contact.hh
  model/solid_mechanics/contact/contact_search_explicit.hh
  model/solid_mechanics/contact/friction_coefficient/unique_constant_fric_coef.hh
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/velocity_weakening_exponential.hh
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/rice_kuwano_modified.hh
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/rice_kuwano.hh
  model/solid_mechanics/contact/friction_coefficient/simplified_dieterich_fric_coef/ruina_slowness_fric_coef.hh
  model/solid_mechanics/contact/friction_coefficient/velocity_weakening_coulomb.hh
  model/solid_mechanics/contact/friction_coefficient/simplified_dieterich_fric_coef.hh
  model/solid_mechanics/contact/friction_coefficient/historic_velocity_fric_coef.hh
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef.hh
  model/solid_mechanics/contact/contact_rigid.hh
  model/solid_mechanics/contact/friction_coefficient.hh
  model/solid_mechanics/contact/contact_search_2d_explicit.hh
  model/solid_mechanics/contact/contact_3d_explicit.hh
  model/solid_mechanics/contact/grid_2d_neighbor_structure.hh
  model/solid_mechanics/contact/regular_grid_neighbor_structure.hh
  model/solid_mechanics/contact/contact_2d_explicit.hh
  model/solid_mechanics/contact_search.hh
  model/solid_mechanics/contact_neighbor_structure.hh

# inline implementation

  model/solid_mechanics/contact/friction_coefficient/unique_constant_fric_coef_inline_impl.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_weakening_coulomb_inline_impl.cc

  model/solid_mechanics/contact/friction_coefficient/unique_constant_fric_coef_inline_impl.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_weakening_coulomb_inline_impl.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/rice_kuwano_modified_inline_impl.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/rice_kuwano_inline_impl.cc
  model/solid_mechanics/contact/friction_coefficient/velocity_dependent_fric_coef/velocity_weakening_exponential_inline_impl.cc
  model/solid_mechanics/contact/friction_coefficient/simplified_dieterich_fric_coef/ruina_slowness_fric_coef_inline_impl.cc
  model/solid_mechanics/contact/friction_coefficient/simplified_dieterich_fric_coef_inline_impl.cc
  model/solid_mechanics/contact/contact_search_explicit_inline_impl.cc
  model/solid_mechanics/contact/contact_search_2d_explicit_inline_impl.cc
  model/solid_mechanics/contact/regular_grid_neighbor_structure_inline_impl.cc
  )


