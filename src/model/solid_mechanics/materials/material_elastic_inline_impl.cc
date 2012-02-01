/**
 * @file   material_elastic_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:57:43 2010
 *
 * @brief  Implementation of the inline functions of the material elastic
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


/* -------------------------------------------------------------------------- */
inline void MaterialElastic::computeStress(Real * F, Real * sigma) {
  Real trace = F[0] + F[4] + F[8]; /// \F_{11} + \F_{22} + \F_{33}

  /// \sigma_{ij} = \lamda * \F_{kk} * \delta_{ij} + 2 * \mu * \F_{ij}
  sigma[0] = lambda * trace + 2*mu*F[0];
  sigma[4] = lambda * trace + 2*mu*F[4];
  sigma[8] = lambda * trace + 2*mu*F[8];

  sigma[1] = sigma[3] =  mu * (F[1] + F[3]);
  sigma[2] = sigma[6] =  mu * (F[2] + F[6]);
  sigma[5] = sigma[7] =  mu * (F[5] + F[7]);
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void  MaterialElastic::computeTangentStiffnessByDim(__attribute__((unused)) akantu::ElementType,
						    akantu::Vector<Real>& tangent_matrix,
						    __attribute__((unused)) akantu::GhostType) {
  AKANTU_DEBUG_IN();

  Real * tangent_val   = tangent_matrix.values;
  UInt offset_tangent  = tangent_matrix.getNbComponent();
  UInt nb_quads        = tangent_matrix.getSize();

  if (nb_quads == 0) return;

  memset(tangent_val, 0, offset_tangent * nb_quads * sizeof(Real));
  for (UInt q = 0; q < nb_quads; ++q, tangent_val += offset_tangent) {
    computeTangentStiffness<dim>(tangent_val);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void  MaterialElastic::computeTangentStiffness(Real * tangent) {

  UInt n = (dim * (dim - 1) / 2 + dim);

  Real Ep = E/((1+nu)*(1-2*nu));
  Real Miiii = Ep * (1-nu);
  Real Miijj = Ep * nu;
  Real Mijij = Ep * (1-2*nu) * .5;

  tangent[0 * n + 0] = Miiii;

  // test of dimension should by optimized out by the compiler due to the template
  if(dim >= 2) {
    tangent[1 * n + 1] = Miiii;
    tangent[0 * n + 1] = Miijj;
    tangent[1 * n + 0] = Miijj;

    tangent[(n - 1) * n + (n - 1)] = Mijij;
  }

  if(dim == 3) {
    tangent[2 * n + 2] = Miiii;
    tangent[0 * n + 2] = Miijj;
    tangent[1 * n + 2] = Miijj;
    tangent[2 * n + 0] = Miijj;
    tangent[2 * n + 1] = Miijj;

    tangent[3 * n + 3] = Mijij;
    tangent[4 * n + 4] = Mijij;
  }
}

/* -------------------------------------------------------------------------- */
inline Real MaterialElastic::getStableTimeStep(Real h, 
					       __attribute__ ((unused)) const Element & element) {
  return (h/getPushWaveSpeed());
}
