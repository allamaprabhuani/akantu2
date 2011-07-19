/**
 * @file   integrator_gauss.cc
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author  Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Feb 10 16:52:07 2011
 *
 * @brief  implementation of gauss integrator class
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
#include "mesh.hh"
#include "integrator_gauss.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
IntegratorGauss::IntegratorGauss(Mesh & mesh,
				 const ID & id,
				 const MemoryID & memory_id) :
  Integrator(mesh, id, memory_id),
  quadrature_points("quadrature_points", id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void IntegratorGauss::checkJacobians(const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = ElementClass<type>::getNbQuadraturePoints();

  UInt * elem_val;
  UInt nb_element;

  elem_val   = mesh->getConnectivity(type,ghost_type).values;
  nb_element = mesh->getConnectivity(type,ghost_type).getSize();

  Real * jacobians_val = jacobians(type, ghost_type).storage();

  for (UInt i = 0; i < nb_element*nb_quadrature_points; ++i,++jacobians_val){
    AKANTU_DEBUG_ASSERT(*jacobians_val >0,
			"Negative jacobian computed,"
			<< " possible problem in the element node order");
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void IntegratorGauss::precomputeJacobiansOnQuadraturePoints(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh->getSpatialDimension();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = ElementClass<type>::getNbQuadraturePoints();

  Real * weights = ElementClass<type>::getGaussIntegrationWeights();

  UInt * elem_val = mesh->getConnectivity(type,ghost_type).values;;
  UInt nb_element = mesh->getConnectivity(type,ghost_type).getSize();

  Vector<Real> & jacobians_tmp = jacobians.alloc(nb_element*nb_quadrature_points,
						 1,
						 type,
						 ghost_type);

  Real * jacobians_val = jacobians_tmp.values;

  Real local_coord[spatial_dimension * nb_nodes_per_element];
  for (UInt elem = 0; elem < nb_element; ++elem) {
    mesh->extractNodalCoordinatesFromElement(local_coord,
					     elem_val+elem*nb_nodes_per_element,
					     nb_nodes_per_element);

    computeJacobianOnQuadPointsByElement<type>(spatial_dimension,
					       local_coord,
					       nb_nodes_per_element,
					       jacobians_val);

    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      *jacobians_val++ *= weights[q];
    }
    //    jacobians_val += nb_quadrature_points;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void IntegratorGauss::integrate(const Vector<Real> & in_f,
				Vector<Real> &intf,
				UInt nb_degre_of_freedom,
				const GhostType & ghost_type,
				const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_element = mesh->getNbElement(type,ghost_type);
  const Vector<Real> & jac_loc = jacobians(type, ghost_type);

  UInt nb_quadrature_points = ElementClass<type>::getNbQuadraturePoints();

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element      = filter_elements->getSize();
    filter_elem_val = filter_elements->values;
  }

  AKANTU_DEBUG_ASSERT(in_f.getSize() == nb_element * nb_quadrature_points,
		      "The vector in_f(" << in_f.getID() << " size " << in_f.getSize()
		      << ") has not the good size (" << nb_element << ").");
  AKANTU_DEBUG_ASSERT(in_f.getNbComponent() == nb_degre_of_freedom ,
		      "The vector in_f(" << in_f.getID()
		      << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(intf.getNbComponent() == nb_degre_of_freedom,
		      "The vector intf(" << intf.getID()
		      << ") has not the good number of component.");


  intf.resize(nb_element);

  Real * in_f_val = in_f.values;
  Real * intf_val = intf.values;
  Real * jac_val  = jac_loc.storage();

  UInt offset_in_f = in_f.getNbComponent()*nb_quadrature_points;
  UInt offset_intf = intf.getNbComponent();

  Real * jac      = jac_val;

  for (UInt el = 0; el < nb_element; ++el) {
    if(filter_elements != NULL) {
      jac      = jac_val  + filter_elem_val[el] * nb_quadrature_points;
    }

    integrate(in_f_val, jac, intf_val, nb_degre_of_freedom, nb_quadrature_points);

    in_f_val += offset_in_f;
    intf_val += offset_intf;
    if(filter_elements == NULL) {
      jac      += nb_quadrature_points;
    }
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template <ElementType type>
Real IntegratorGauss::integrate(const Vector<Real> & in_f,
				const GhostType & ghost_type,
				const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_element = mesh->getNbElement(type, ghost_type);
  const Vector<Real> & jac_loc = jacobians(type, ghost_type);

  UInt nb_quadrature_points = ElementClass<type>::getNbQuadraturePoints();

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element      = filter_elements->getSize();
    filter_elem_val = filter_elements->values;
  }

  AKANTU_DEBUG_ASSERT(in_f.getSize() == nb_element * nb_quadrature_points,
		      "The vector in_f(" << in_f.getID()
		      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(in_f.getNbComponent() == 1,
		      "The vector in_f(" << in_f.getID()
		      << ") has not the good number of component.");

  Real intf = 0.;
  Real * in_f_val  = in_f.values;
  Real * jac_val   = jac_loc.storage();
  UInt offset_in_f = in_f.getNbComponent() * nb_quadrature_points;
  Real * jac       = jac_val;

  for (UInt el = 0; el < nb_element; ++el) {
    if(filter_elements != NULL) {
      jac = jac_val  + filter_elem_val[el] * nb_quadrature_points;
    }

    Real el_intf = 0;
    integrate(in_f_val, jac, &el_intf, 1, nb_quadrature_points);
    intf += el_intf;

    in_f_val += offset_in_f;
    if(filter_elements == NULL) {
      jac += nb_quadrature_points;
    }
  }

  AKANTU_DEBUG_OUT();
  return intf;
}



/* -------------------------------------------------------------------------- */
template <>
void IntegratorGauss::precomputeJacobiansOnQuadraturePoints<_bernoulli_beam_2>(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh->getSpatialDimension();

  UInt nb_quadrature_points = ElementClass<_bernoulli_beam_2>::getNbQuadraturePoints();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(_bernoulli_beam_2);

  UInt * elem_val;
  UInt nb_element;

  elem_val   = mesh->getConnectivity(_bernoulli_beam_2, ghost_type).values;
  nb_element = mesh->getConnectivity(_bernoulli_beam_2, ghost_type).getSize();

  Vector<Real> & jacobians_tmp = jacobians.alloc(nb_element*nb_quadrature_points,
						 1, _bernoulli_beam_2, ghost_type);

  Real local_coord[spatial_dimension * nb_nodes_per_element];
  Real * jacobians_val = jacobians_tmp.values;
  for (UInt elem = 0; elem < nb_element; ++elem) {
    mesh->extractNodalCoordinatesFromElement(local_coord,
					     elem_val+elem*nb_nodes_per_element,
					     nb_nodes_per_element);

    ElementClass<_bernoulli_beam_2>::computeJacobian(local_coord,
						     nb_quadrature_points,
						     spatial_dimension,
						     jacobians_val);
    jacobians_val += nb_quadrature_points;
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/* template instanciation */
/* -------------------------------------------------------------------------- */

#define INSTANCIATE_TEMPLATE_CLASS(type)				\
  template void IntegratorGauss::					\
  precomputeJacobiansOnQuadraturePoints<type>(const GhostType & ghost_type); \
  template void IntegratorGauss::					\
  checkJacobians<type>(const GhostType & ghost_type) const;		\
  template void IntegratorGauss::					\
  integrate<type>(const Vector<Real> & in_f,				\
		  Vector<Real> &intf,					\
		  UInt nb_degre_of_freedom,				\
		  const GhostType & ghost_type,				\
		  const Vector<UInt> * filter_elements) const;		\
  template Real IntegratorGauss::					\
  integrate<type>(const Vector<Real> & in_f,				\
		  const GhostType & ghost_type,				\
		  const Vector<UInt> * filter_elements) const;


AKANTU_BOOST_ELEMENT_LIST(INSTANCIATE_TEMPLATE_CLASS)
#undef INSTANCIATE_TEMPLATE_CLASS


__END_AKANTU__
