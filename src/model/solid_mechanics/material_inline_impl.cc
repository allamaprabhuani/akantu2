/**
 * @file   material_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:57:43 2010
 *
 * @brief  Implementation of the inline functions of the class material
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
inline UInt Material::addElement(const ElementType & type,
				 UInt element,
				 const GhostType & ghost_type) {
  element_filter(type, ghost_type).push_back(element);
  return element_filter(type, ghost_type).getSize()-1;
}

/* -------------------------------------------------------------------------- */
inline UInt Material::getTangentStiffnessVoigtSize(UInt dim) const {
  return (dim * (dim - 1) / 2 + dim);
}

/* -------------------------------------------------------------------------- */
inline void Material::computePotentialEnergy(Real * F, Real * sigma, Real * epot) {
  *epot = 0.;
  for (UInt i = 0, t = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j, ++t)
      (*epot) += sigma[t] * F[t];

  *epot *= .5;
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void Material::transferBMatrixToSymVoigtBMatrix(Real * B, Real * Bvoigt, UInt nb_nodes_per_element) const {
  UInt size = getTangentStiffnessVoigtSize(dim) * nb_nodes_per_element * dim;
  memset(Bvoigt, 0, size * sizeof(Real));

  for (UInt i = 0; i < dim; ++i) {
    Real * Bvoigt_tmp = Bvoigt + i * (dim * nb_nodes_per_element + 1);
    Real * Btmp = B + i;
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      *Bvoigt_tmp = *Btmp;
      Btmp += dim;
      Bvoigt_tmp += dim;
    }
  }

  if(dim == 2) {
    ///in 2D, fill the @f$ [\frac{\partial N_i}{\partial x}, \frac{\partial N_i}{\partial y}]@f$ row
    Real * Bvoigt_tmp = Bvoigt + dim * nb_nodes_per_element * 2;
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      Bvoigt_tmp[1] = B[n * dim + 0];
      Bvoigt_tmp[0] = B[n * dim + 1];
      Bvoigt_tmp += dim;
    }
  }


  if(dim == 3) {
    UInt Bvoigt_wcol = dim * nb_nodes_per_element;
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      Real dndx = B[n * dim + 0];
      Real dndy = B[n * dim + 1];
      Real dndz = B[n * dim + 2];

      UInt Bni_off = n * dim;

      ///in 3D, fill the @f$ [0, \frac{\partial N_i}{\partial y}, \frac{N_i}{\partial z}]@f$ row
      Bvoigt[3 * Bvoigt_wcol + Bni_off + 1] = dndz;
      Bvoigt[3 * Bvoigt_wcol + Bni_off + 2] = dndy;

      ///in 3D, fill the @f$ [\frac{\partial N_i}{\partial x}, 0, \frac{N_i}{\partial z}]@f$ row
      Bvoigt[4 * Bvoigt_wcol + Bni_off + 0] = dndz;
      Bvoigt[4 * Bvoigt_wcol + Bni_off + 2] = dndx;

      ///in 3D, fill the @f$ [\frac{\partial N_i}{\partial x}, \frac{N_i}{\partial y}, 0]@f$ row
      Bvoigt[5 * Bvoigt_wcol + Bni_off + 0] = dndy;
      Bvoigt[5 * Bvoigt_wcol + Bni_off + 1] = dndx;
    }
  }
}

