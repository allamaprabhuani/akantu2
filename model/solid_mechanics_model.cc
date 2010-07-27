/**
 * @file   solid_mechanics_model.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 22 14:35:38 2010
 *
 * @brief  Implementation of the SolidMechanicsModel class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
SolidMechanicsModel::SolidMechanicsModel(UInt spatial_dimension,
					 const ModelID & id,
					 const MemoryID & memory_id) : 
  Model(spatial_dimension, id, memory_id) {
  AKANTU_DEBUG_IN();
  
  for(UInt t = _not_defined; t < _max_element_type; ++t) {
    this->displacement    [t] = NULL;
    this->mass            [t] = NULL;
    this->velocity        [t] = NULL;
    this->acceleration    [t] = NULL;
    this->force           [t] = NULL;
    this->residual        [t] = NULL;
    this->stress          [t] = NULL;
    this->strain          [t] = NULL;
    this->boundary        [t] = NULL;
    this->element_material[t] = NULL;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolidMechanicsModel::SolidMechanicsModel(UInt spatial_dimension,
					 Mesh & mesh,
					 const ModelID & id,
					 const MemoryID & memory_id) :
  Model(spatial_dimension, mesh, id, memory_id) {
  AKANTU_DEBUG_IN();

  for(UInt t = _not_defined; t < _max_element_type; ++t) {
    this->displacement    [t] = NULL;
    this->mass            [t] = NULL;
    this->velocity        [t] = NULL;
    this->acceleration    [t] = NULL;
    this->force           [t] = NULL;
    this->residual        [t] = NULL;
    this->stress          [t] = NULL;
    this->strain          [t] = NULL;
    this->boundary        [t] = NULL;
    this->element_material[t] = NULL;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolidMechanicsModel::~SolidMechanicsModel() {
  const Mesh::ConnectivityTypeList & type_list = fem->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    AKANTU_DEBUG(dblAccessory, "Deleting displacements vector of type " << *it);
    dealloc(displacement[*it]->getID());
    displacement[*it] = NULL;

    AKANTU_DEBUG(dblAccessory, "Deleting mass vector of type " << *it);
    dealloc(mass[*it]->getID());
    mass[*it] = NULL;

    AKANTU_DEBUG(dblAccessory, "Deleting velocity vector of type " << *it);
    dealloc(velocity[*it]->getID());
    velocity[*it] = NULL;

    AKANTU_DEBUG(dblAccessory, "Deleting acceleration vector of type " << *it);
    dealloc(acceleration[*it]->getID());
    acceleration[*it] = NULL;

    AKANTU_DEBUG(dblAccessory, "Deleting force vector of type " << *it);
    dealloc(force[*it]->getID());
    force[*it] = NULL;

    AKANTU_DEBUG(dblAccessory, "Deleting residual vector of type " << *it);
    dealloc(residual[*it]->getID());
    residual[*it] = NULL;

    AKANTU_DEBUG(dblAccessory, "Deleting stress derivatives vector of type " << *it);
    dealloc(stress[*it]->getID());
    stress[*it] = NULL;

    AKANTU_DEBUG(dblAccessory, "Deleting strain vector of type " << *it);
    dealloc(strain[*it]->getID());
    strain[*it] = NULL;

    AKANTU_DEBUG(dblAccessory, "Deleting boundary vector of type " << *it);
    dealloc(boundary[*it]->getID());
    boundary[*it] = NULL;

    AKANTU_DEBUG(dblAccessory, "Deleting element_material vector of type " << *it);
    dealloc(element_material[*it]->getID());
    element_material[*it] = NULL;
  }

}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initModel() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = fem->getNbNodes();

  const Mesh::ConnectivityTypeList & type_list = fem->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    UInt nb_quadrature_points = fem->getNbQuadraturePoints(*it);
    UInt nb_element           = fem->getNbElement(*it);

    std::stringstream sstr_disp; sstr_disp << id << ":displacement:" << *it;
    std::stringstream sstr_mass; sstr_mass << id << ":mass:" << *it;
    std::stringstream sstr_velo; sstr_velo << id << ":velocity:" << *it;
    std::stringstream sstr_acce; sstr_acce << id << ":acceleration:" << *it;
    std::stringstream sstr_forc; sstr_forc << id << ":force:" << *it;
    std::stringstream sstr_resi; sstr_resi << id << ":residual:" << *it;
    std::stringstream sstr_stre; sstr_stre << id << ":stress:" << *it;
    std::stringstream sstr_stra; sstr_stra << id << ":strain:" << *it;
    std::stringstream sstr_boun; sstr_boun << id << ":boundary:" << *it;
    std::stringstream sstr_elma; sstr_elma << id << ":element_material:" << *it;
    displacement    [*it] = &(alloc<Real>(sstr_disp.str(), nb_nodes,
					  spatial_dimension, NAN));
    mass            [*it] = &(alloc<Real>(sstr_mass.str(), nb_nodes,
					  1)); // \todo see if it must not be spatial_dimension
    velocity        [*it] = &(alloc<Real>(sstr_velo.str(), nb_nodes,
					  spatial_dimension, NAN));
    acceleration    [*it] = &(alloc<Real>(sstr_acce.str(), nb_nodes,
					  spatial_dimension, NAN));
    force           [*it] = &(alloc<Real>(sstr_forc.str(), nb_nodes,
					  spatial_dimension, NAN));
    residual        [*it] = &(alloc<Real>(sstr_resi.str(), nb_nodes,
					  spatial_dimension, NAN));
    stress          [*it] = &(alloc<Real>(sstr_stre.str(), nb_element,
					  nb_quadrature_points * spatial_dimension * spatial_dimension,
					  NAN));
    strain          [*it] = &(alloc<Real>(sstr_stra.str(), nb_element,
					  nb_quadrature_points * spatial_dimension * spatial_dimension,
					  NAN));
    boundary        [*it] = &(alloc<Int>(sstr_boun.str(), nb_nodes,
					 spatial_dimension, 0));
    element_material[*it] = &(alloc<Int>(sstr_elma.str(), nb_element, 1, 0));
  }

  std::stringstream sstr; sstr << id << ":0:elastic" << *it;
  materials.push_back(Material<_elastic>(*this, sstr.str()));

  fem->initShapeFunctions();

  AKANTU_DEBUG_OUT();
};

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMass() {
  AKANTU_DEBUG_IN();

  MaterialBase * mat_val = &(materials.at(0));

  const Mesh::ConnectivityTypeList & type_list = fem->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(fem->getSpatialDimension(*it) != spatial_dimension) continue;

    UInt nb_nodes_per_element = fem->getNbNodesPerElement(*it);
    UInt nb_element           = fem->getNbElement(*it);
    
    const Vector<Real> & shapes = fem->getShapes(*it);
    Vector<Real> * rho_phi_i = new Vector<Real>(nb_element, nb_nodes_per_element);

    Int * elem_mat_val = element_material[*it]->values;
    Real * rho_phi_i_val = rho_phi_i->values;
    Real * shapes_val = shapes.values;

    for (UInt el = 0; el < nb_element; ++el) {
      UInt rho = mat_val[elem_mat_val[el]].getRho();
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	*rho_phi_i_val++ = rho * *shapes_val++;
      }
    }

    Vector<Real> * int_rho_phi_i = new Vector<Real>(nb_element, nb_nodes_per_element);
    fem->integrate(*rho_phi_i, *int_rho_phi_i, nb_nodes_per_element, *it);
    delete rho_phi_i;

    fem->assembleVector(*int_rho_phi_i, *mass[*it], 1, *it);
    delete int_rho_phi_i;
  }

  
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::readMaterials(const std::string & filename) {

}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  
  stream << space << "Solid Mechanics Model [" << std::endl;
  stream << space << " + id                : " << id << std::endl;
  stream << space << " + spatial dimension : " << spatial_dimension << std::endl;

  stream << space << " + fem [" << std::endl;
  fem->printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;
  stream << space << " + connectivity type information [" << std::endl;
  const Mesh::ConnectivityTypeList & type_list = fem->getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    stream << space << AKANTU_INDENT << AKANTU_INDENT << " + " << *it <<" [" << std::endl;
    displacement    [*it]->printself(stream, indent + 3);
    mass            [*it]->printself(stream, indent + 3);
    velocity        [*it]->printself(stream, indent + 3);
    acceleration    [*it]->printself(stream, indent + 3);
    force           [*it]->printself(stream, indent + 3);
    residual        [*it]->printself(stream, indent + 3);
    stress          [*it]->printself(stream, indent + 3);
    strain          [*it]->printself(stream, indent + 3);
    boundary        [*it]->printself(stream, indent + 3);
    element_material[*it]->printself(stream, indent + 3);
    stream << space << AKANTU_INDENT << AKANTU_INDENT << "]" << std::endl;
  }
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
