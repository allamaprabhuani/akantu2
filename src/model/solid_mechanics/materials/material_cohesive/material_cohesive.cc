/**
 * @file   material_cohesive.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Tue Feb  7 18:24:52 2012
 *
 * @brief  Specialization of the material class for cohesive elements
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
#include "material_cohesive.hh"
#include "solid_mechanics_model.hh"
#include "sparse_matrix.hh"
#include "dof_synchronizer.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialCohesive::MaterialCohesive(Model & model, const ID & id) :
  Material(model,id),
  reversible_energy("reversible_energy", id),
  total_energy("total_energy", id),
  tractions_old("tractions (old)",id),
  opening_old("opening (old)",id),
  tractions("tractions",id),
  opening("opening",id) {

  AKANTU_DEBUG_IN();

  initInternalVector(reversible_energy, 1, _ek_cohesive);
  initInternalVector(total_energy, 1, _ek_cohesive);
  initInternalVector(tractions_old, spatial_dimension, _ek_cohesive);
  initInternalVector(tractions, spatial_dimension, _ek_cohesive);
  initInternalVector(opening_old, spatial_dimension, _ek_cohesive);
  initInternalVector(opening, spatial_dimension, _ek_cohesive);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MaterialCohesive::~MaterialCohesive() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool MaterialCohesive::setParam(const std::string & key, const std::string & value,
				const ID & id) {
  return Material::setParam(key,value,id);
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::initMaterial() {
  AKANTU_DEBUG_IN();

  Material::initMaterial();

  fem_cohesive = &(model->getFEMClass<MyFEMCohesiveType>("CohesiveFEM"));

  Mesh & mesh = fem_cohesive->getMesh();

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    Mesh::type_iterator it = mesh.firstType(spatial_dimension, gt, _ek_cohesive);
    Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, gt, _ek_cohesive);

    for(; it != last_type; ++it) {
      ElementType type = *it;
      element_filter.alloc(0, 1, type);
    }
  }

  resizeCohesiveVectors();

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void MaterialCohesive::resizeCohesiveVectors() {

  resizeInternalVector(reversible_energy, _ek_cohesive);
  resizeInternalVector(total_energy, _ek_cohesive);
  resizeInternalVector(tractions_old, _ek_cohesive);
  resizeInternalVector(tractions, _ek_cohesive);
  resizeInternalVector(opening_old, _ek_cohesive);
  resizeInternalVector(opening, _ek_cohesive);

}

/* -------------------------------------------------------------------------- */
/**
 * Compute  the  residual  by  assembling  @f$\int_{e}  t_e  N_e dS @f$
 *
 * @param[in] displacements nodes displacements
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void MaterialCohesive::updateResidual(__attribute__((unused)) Vector<Real> & displacement,
				      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// compute traction
  computeTraction(ghost_type);

  /// update and assemble residual
  assembleResidual(ghost_type);

  /// compute energies
  computeEnergies();

  /// update old values
  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type, _ek_cohesive);
  for(; it != last_type; ++it) {
    tractions_old(*it, ghost_type).copy(tractions(*it, ghost_type));
    opening_old(*it, ghost_type).copy(opening(*it, ghost_type));
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void MaterialCohesive::assembleResidual(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model->getSpatialDimension();

  Vector<Real> & residual = const_cast<Vector<Real> &>(model->getResidual());

  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type, _ek_cohesive);
  for(; it != last_type; ++it) {
    const Vector<Real> & shapes = fem_cohesive->getShapes(*it, ghost_type);

    Vector<UInt> & elem_filter = element_filter(*it, ghost_type);
    Vector<Real> & traction = tractions(*it, ghost_type);

    UInt size_of_shapes       = shapes.getNbComponent();
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
    UInt nb_quadrature_points = fem_cohesive->getNbQuadraturePoints(*it, ghost_type);

    UInt nb_element = elem_filter.getSize();

    /// compute @f$t_i N_a@f$

    Real * shapes_val       = shapes.storage();
    UInt * elem_filter_val  = elem_filter.storage();

    Vector<Real> * shapes_filtered =
      new Vector<Real>(nb_element*nb_quadrature_points, size_of_shapes, "filtered shapes");
    Real * shapes_filtered_val = shapes_filtered->values;

    for (UInt el = 0; el < nb_element; ++el) {
      shapes_val = shapes.storage() + elem_filter_val[el] *
	size_of_shapes * nb_quadrature_points;
      memcpy(shapes_filtered_val, shapes_val,
	     size_of_shapes * nb_quadrature_points * sizeof(Real));
      shapes_filtered_val += size_of_shapes * nb_quadrature_points;
    }

    shapes_filtered_val = shapes_filtered->values;

    // multiply traction by shapes

    Vector<Real> traction_cpy(traction);
    traction_cpy.extendComponentsInterlaced(size_of_shapes, spatial_dimension);

    Real * traction_cpy_val = traction_cpy.storage();

    for (UInt el = 0; el < nb_element; ++el) {
      for (UInt q = 0; q < nb_quadrature_points; ++q) {
	for (UInt n = 0; n < size_of_shapes; ++n,++shapes_filtered_val) {
	  for (UInt i = 0; i < spatial_dimension; ++i) {
	    *traction_cpy_val++ *= *shapes_filtered_val;
	  }
	}
      }
    }

    delete shapes_filtered;

    /**
     * compute @f$\int t \cdot N\, dS@f$ by  @f$ \sum_q \mathbf{N}^t
     * \mathbf{t}_q \overline w_q J_q@f$
     */
    Vector<Real> int_t_N(nb_element,spatial_dimension*size_of_shapes,
			 "int_t_N");

    fem_cohesive->integrate(traction_cpy, int_t_N,
			    spatial_dimension*size_of_shapes,
			    *it, ghost_type,
			    &elem_filter);

    int_t_N.extendComponentsInterlaced(2,spatial_dimension);

    Real * int_t_N_val = int_t_N.storage();
    for (UInt el = 0; el < nb_element; ++el) {
      for (UInt n = 0; n < size_of_shapes*spatial_dimension; ++n) {
	int_t_N_val[n] *= -1.;
      }
      int_t_N_val += nb_nodes_per_element*spatial_dimension;
    }


    /// assemble
    model->getFEMBoundary().assembleVector(int_t_N, residual,
					   model->getDOFSynchronizer().getLocalDOFEquationNumbers(),
					   residual.getNbComponent(),
					   *it, ghost_type, &elem_filter, -1);

  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute traction from displacements
 *
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void MaterialCohesive::computeTraction(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model->getSpatialDimension();

  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type, _ek_cohesive);

  for(; it != last_type; ++it) {
    Vector<UInt> & elem_filter = element_filter(*it, ghost_type);

    UInt nb_element = elem_filter.getSize();
    UInt nb_quadrature_points = nb_element*fem_cohesive->getNbQuadraturePoints(*it, ghost_type);

    Vector<Real> normal(nb_quadrature_points, spatial_dimension,
			"normal");

    /// compute normals @f$\mathbf{n}@f$
    computeNormal(model->getCurrentPosition(), normal, *it, ghost_type);

    /// compute openings @f$\mathbf{\delta}@f$
    computeOpening(model->getDisplacement(), opening(*it, ghost_type), *it, ghost_type);

    /// compute traction @f$\mathbf{t}@f$
    computeTraction(normal, *it, ghost_type);

  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void MaterialCohesive::computeNormal(const Vector<Real> & position,
				     Vector<Real> & normal,
				     ElementType type,
				     GhostType ghost_type) {
  AKANTU_DEBUG_IN();

#define COMPUTE_NORMAL(type)						\
  fem_cohesive->getShapeFunctions().					\
    computeNormalsOnControlPoints<type, CohesiveReduceFunctionMean>(position, \
								    normal, \
								    ghost_type,	\
								    &(element_filter(type, ghost_type)))

  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_NORMAL);
#undef COMPUTE_NORMAL

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void MaterialCohesive::computeOpening(const Vector<Real> & displacement,
				      Vector<Real> & opening,
				      ElementType type,
				      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

#define COMPUTE_OPENING(type)						\
  fem_cohesive->getShapeFunctions().					\
    interpolateOnControlPoints<type, CohesiveReduceFunctionOpening>(displacement, \
								    opening, \
								    spatial_dimension, \
								    ghost_type, \
								    &(element_filter(type, ghost_type)))

  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_OPENING);
#undef COMPUTE_OPENING

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::computeEnergies() {
  AKANTU_DEBUG_IN();

  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _not_ghost, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _not_ghost, _ek_cohesive);

  for(; it != last_type; ++it) {
    Vector<Real>::iterator<Real> erev =
      reversible_energy(*it, _not_ghost).begin();
    Vector<Real>::iterator<Real> etot =
      total_energy(*it, _not_ghost).begin();
    Vector<Real>::iterator<types::RVector> traction_it =
      tractions(*it, _not_ghost).begin(spatial_dimension);
    Vector<Real>::iterator<types::RVector> traction_old_it =
      tractions_old(*it, _not_ghost).begin(spatial_dimension);
    Vector<Real>::iterator<types::RVector> opening_it =
      opening(*it, _not_ghost).begin(spatial_dimension);
    Vector<Real>::iterator<types::RVector> opening_old_it =
      opening_old(*it, _not_ghost).begin(spatial_dimension);

    Vector<Real>::iterator<types::RVector>traction_end =
      tractions(*it, _not_ghost).end(spatial_dimension);

    /// loop on each quadrature point
    for (; traction_it != traction_end;
	 ++traction_it, ++traction_old_it,
	   ++opening_it, ++opening_old_it,
	   ++erev, ++etot) {

      /// trapezoidal integration
      types::RVector b(spatial_dimension);
      b  = *opening_it;
      b -= *opening_old_it;

      types::RVector h(spatial_dimension);
      h  = *traction_old_it;
      h += *traction_it;

      *etot += .5 * b.dot(h);
      *erev  = .5 * traction_it->dot(*opening_it);
    }


  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getReversibleEnergy() {
  AKANTU_DEBUG_IN();
  Real erev = 0.;

  /// integrate the dissipated energy for each type of elements
  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _not_ghost, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _not_ghost, _ek_cohesive);

  for(; it != last_type; ++it) {
    erev += fem_cohesive->integrate(reversible_energy(*it, _not_ghost), *it,
				    _not_ghost, &element_filter(*it, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return erev;
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getDissipatedEnergy() {
  AKANTU_DEBUG_IN();
  Real edis = 0.;

  /// integrate the dissipated energy for each type of elements
  Mesh & mesh = fem_cohesive->getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, _not_ghost, _ek_cohesive);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, _not_ghost, _ek_cohesive);

  for(; it != last_type; ++it) {
    Vector<Real> dissipated_energy(total_energy(*it, _not_ghost));
    dissipated_energy -= reversible_energy(*it, _not_ghost);
    edis += fem_cohesive->integrate(dissipated_energy, *it,
				    _not_ghost, &element_filter(*it, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return edis;
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Material Cohesive [" << std::endl;
  Material::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}


__END_AKANTU__
