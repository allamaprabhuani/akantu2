#include "material_bulk_viscosity.hh"

namespace akantu {

template <Int dim>
using MaterialElasticBulkViscosity =
    class MaterialBulkViscosity<dim, MaterialElastic<dim>>;

static bool material_is_alocated_elastic_bulk_viscosity =
    instantiateMaterial<MaterialElasticBulkViscosity>("elastic_bulk_viscosity");

} // namespace akantu
