/**
 * @file   weight_function.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Mar 23 15:55:58 2012
 *
 * @brief  implementation of the weight function classes
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

#include "weight_function.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* Stress based weight function                                               */
/* -------------------------------------------------------------------------- */
StressBasedWeightFunction::StressBasedWeightFunction(Real radius, Real ft,
						     UInt spatial_dimension,
						     const Material & material) :
  BaseWeightFunction(radius), material(&material), ft(ft),
  spatial_dimension(spatial_dimension),
  stress_diag("stress_diag", material.getID()), selected_stress_diag(NULL),
  stress_base("stress_base", material.getID()), selected_stress_base(NULL),
  characteristic_size("lc", material.getID()),  selected_characteristic_size(NULL)
{
  material.initInternalVector(stress_diag, spatial_dimension);
  material.resizeInternalVector(stress_diag);
  material.initInternalVector(stress_base, spatial_dimension * spatial_dimension);
  material.resizeInternalVector(stress_base);
  material.initInternalVector(characteristic_size, 1);
  material.resizeInternalVector(characteristic_size);

  const Mesh & mesh = material.getModel().getFEM().getMesh();
  for (UInt g = _not_ghost; g < _casper; ++g) {
    GhostType gt = GhostType(g);
    Mesh::type_iterator it = mesh.firstType(spatial_dimension, gt);
    Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, gt);
    for(; it != last_type; ++it) {
      UInt nb_quadrature_points =
	material.getModel().getFEM().getNbQuadraturePoints(*it, gt);
      const Vector<UInt> & element_filter = material.getElementFilter(*it, gt);
      UInt nb_element = element_filter.getSize();

      Vector<Real> ones(nb_element*nb_quadrature_points, 1, 1.);
      Vector<Real> & lc = characteristic_size(*it, gt);
      material.getModel().getFEM().integrateOnQuadraturePoints(ones,
							       lc,
							       1,
							       *it,
							       gt,
							       &element_filter);

      for (UInt q = 0;  q < nb_quadrature_points * nb_element; q++) {
	lc(q) = pow(lc(q), 1./ Real(spatial_dimension));
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void StressBasedWeightFunction::updatePrincipalStress(GhostType ghost_type) {
  const Mesh & mesh = material->getModel().getFEM().getMesh();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    Vector<Real>::const_iterator<types::Matrix> sigma =
      material->getStress(*it, ghost_type).begin(spatial_dimension, spatial_dimension);
    Vector<Real>::iterator<types::RVector> eigenvalues =
      stress_diag(*it, ghost_type).begin(spatial_dimension);
    Vector<Real>::iterator<types::RVector> eigenvalues_end =
      stress_diag(*it, ghost_type).end(spatial_dimension);
    Vector<Real>::iterator<types::Matrix> eigenvector =
      stress_base(*it, ghost_type).begin(spatial_dimension, spatial_dimension);

#ifndef __trick__
    Vector<Real>::iterator<Real> cl = characteristic_size(*it, ghost_type).begin();
#endif
    for(;eigenvalues != eigenvalues_end; ++eigenvalues, ++eigenvector, ++cl) {
      sigma->eig(*eigenvalues, *eigenvector);
      *eigenvalues /= ft;

#ifndef __trick__
      // specify a lower bound for principal stress based on the size of the element
      for (UInt i = 0; i < spatial_dimension; ++i) {
        (*eigenvalues)(i) = std::max(*cl / R, (*eigenvalues)(i));
      }
#endif
    }
  }
}


__END_AKANTU__
