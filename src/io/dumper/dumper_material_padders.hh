/**
 * Copyright (©) 2014-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef AKANTU_DUMPER_MATERIAL_PADDERS_HH_
#define AKANTU_DUMPER_MATERIAL_PADDERS_HH_
/* -------------------------------------------------------------------------- */
#include "dumper_padding_helper.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
namespace dumpers {
  /* ------------------------------------------------------------------------ */
  class MaterialFunctor {
    /* ---------------------------------------------------------------------- */
    /* Constructors/Destructors                                               */
    /* ---------------------------------------------------------------------- */
  public:
    MaterialFunctor(const SolidMechanicsModel & model)
        : model(model), material_index(model.getMaterialByElement()),
          nb_data_per_element("nb_data_per_element", model.getID()),
          spatial_dimension(model.getSpatialDimension()) {}

    /* ---------------------------------------------------------------------- */
    /* Methods                                                                */
    /* ---------------------------------------------------------------------- */
    /// return the material from the global element index
    const Material & getMaterialFromGlobalIndex(Element global_index) {
      auto index = global_index.element;
      auto material_id = material_index(global_index.type)(index);
      const Material & material = model.getMaterial(material_id);
      return material;
    }

    /// return the type of the element from global index
    ElementType
    getElementTypeFromGlobalIndex( // NOLINT(readability-convert-member-functions-to-static)
        Element global_index) {
      return global_index.type;
    }

  protected:
    /* ---------------------------------------------------------------------- */
    /* Class Members                                                          */
    /* ---------------------------------------------------------------------- */

    /// all material padders probably need access to solid mechanics model
    const SolidMechanicsModel & model;

    /// they also need an access to the map from global ids to material id and
    /// local ids
    const ElementTypeMapArray<Idx> & material_index;

    /// the number of data per element
    const ElementTypeMapArray<Idx> nb_data_per_element;

    Int spatial_dimension;
  };

  /* ------------------------------------------------------------------------ */
  template <class T, class R>
  class MaterialPadder : public MaterialFunctor,
                         public PadderGeneric<Vector<T>, R> {
  public:
    MaterialPadder(const SolidMechanicsModel & model)
        : MaterialFunctor(model) {}
  };

  /* ------------------------------------------------------------------------ */

  template <Int spatial_dimension>
  class StressPadder : public MaterialPadder<Real, Matrix<Real>> {

  public:
    StressPadder(const SolidMechanicsModel & model)
        : MaterialPadder<Real, Matrix<Real>>(model) {
      this->setPadding(3, 3);
    }

    inline Matrix<Real> func(const Vector<Real> & in,
                             Element global_element_id) override {
      auto nrows = spatial_dimension;
      auto ncols = in.size() / nrows;
      auto nb_data = in.size() / (nrows * nrows);

      Matrix<Real> stress = this->pad(in, nrows, ncols, nb_data);
      const Material & material =
          this->getMaterialFromGlobalIndex(global_element_id);
      bool plane_strain = true;
      if (spatial_dimension == 2) {
        plane_strain = !((bool)material.getParam("Plane_Stress"));
      }

      if (plane_strain) {
        Real nu = material.getParam("nu");
        for (Int d = 0; d < nb_data; ++d) {
          stress(2, 2 + 3 * d) =
              nu * (stress(0, 0 + 3 * d) + stress(1, 1 + 3 * d));
        }
      }
      return stress;
    }

    Int getDim() override { return 9; }
    Int getNbComponent(Int /*old_nb_comp*/) override { return this->getDim(); }
  };

  /* ------------------------------------------------------------------------ */
  template <Int spatial_dimension>
  class StrainPadder : public MaterialFunctor,
                       public PadderGeneric<Matrix<Real>, Matrix<Real>> {
  public:
    StrainPadder(const SolidMechanicsModel & model) : MaterialFunctor(model) {
      this->setPadding(3, 3);
    }

    inline Matrix<Real> func(const Matrix<Real> & in,
                             Element global_element_id) override {
      auto nrows = spatial_dimension;
      auto nb_data = in.size() / (nrows * nrows);

      Matrix<Real> strain = this->pad(in, nb_data);
      const Material & material =
          this->getMaterialFromGlobalIndex(global_element_id);
      bool plane_stress = material.getParam("Plane_Stress");

      if (plane_stress) {
        Real nu = material.getParam("nu");
        for (Int d = 0; d < nb_data; ++d) {
          strain(2, 2 + 3 * d) =
              nu / (nu - 1) * (strain(0, 0 + 3 * d) + strain(1, 1 + 3 * d));
        }
      }

      return strain;
    }

    Int getDim() override { return 9; }
    Int getNbComponent(Int /*old_nb_comp*/) override { return this->getDim(); }
  };

  /* ------------------------------------------------------------------------ */
  template <bool green_strain>
  class ComputeStrain : public MaterialFunctor,
                        public ComputeFunctor<Vector<Real>, Matrix<Real>> {
  public:
    ComputeStrain(const SolidMechanicsModel & model) : MaterialFunctor(model) {}

    inline Matrix<Real> func(const Vector<Real> & in,
                             Element /*global_element_id*/) override {
      auto nrows = spatial_dimension;
      auto ncols = in.size() / nrows;
      auto nb_data = in.size() / (nrows * nrows);

      Matrix<Real> ret_all_strain(nrows, ncols);
      Tensor3Proxy<const Real> all_grad_u(in.data(), nrows, nrows, nb_data);
      Tensor3Proxy<Real> all_strain(ret_all_strain.data(), nrows, nrows,
                                    nb_data);

      for (Int d = 0; d < nb_data; ++d) {
        auto && grad_u = all_grad_u(d);
        auto && strain = all_strain(d);

        if (spatial_dimension == 2) {
          if (green_strain) {
            strain = Material::gradUToE<2>(grad_u);
          } else {
            strain = Material::gradUToEpsilon<2>(grad_u);
          }
        } else if (spatial_dimension == 3) {
          if (green_strain) {
            strain = Material::gradUToE<3>(grad_u);
          } else {
            strain = Material::gradUToEpsilon<3>(grad_u);
          }
        }
      }

      return ret_all_strain;
    }

    Int getDim() override { return spatial_dimension * spatial_dimension; }
    Int getNbComponent(Int /*old_nb_comp*/) override { return this->getDim(); }
  };

  /* ------------------------------------------------------------------------ */
  template <bool green_strain>
  class ComputePrincipalStrain
      : public MaterialFunctor,
        public ComputeFunctor<Vector<Real>, Matrix<Real>> {
  public:
    ComputePrincipalStrain(const SolidMechanicsModel & model)
        : MaterialFunctor(model) {}

    inline Matrix<Real> func(const Vector<Real> & in,
                             Element /*global_element_id*/) override {
      auto nrows = spatial_dimension;
      auto nb_data = in.size() / (nrows * nrows);

      Matrix<Real> ret_all_strain(nrows, nb_data);
      Tensor3Proxy<const Real> all_grad_u(in.data(), nrows, nrows, nb_data);
      Matrix<Real> strain(nrows, nrows);

      for (Int d = 0; d < nb_data; ++d) {
        Matrix<Real> grad_u = all_grad_u(d);

        tuple_dispatch<AllSpatialDimensions>(
            [&grad_u, &strain](auto && dim_t) {
              constexpr auto dim = aka::decay_v<decltype(dim_t)>;
              if (green_strain) {
                Material::gradUToE<dim>(grad_u, strain);
              } else {
                strain = Material::gradUToEpsilon<dim>(grad_u);
              }
            },
            spatial_dimension);

        auto && principal_strain = ret_all_strain(d);
        strain.eig(principal_strain);
      }

      return ret_all_strain;
    }

    Int getDim() override { return spatial_dimension; }
    Int getNbComponent(Int /*old_nb_comp*/) override { return this->getDim(); }
  };

  /* ------------------------------------------------------------------------ */
  class ComputeVonMisesStress
      : public MaterialFunctor,
        public ComputeFunctor<Vector<Real>, Vector<Real>> {
  public:
    ComputeVonMisesStress(const SolidMechanicsModel & model)
        : MaterialFunctor(model) {}

    inline Vector<Real> func(const Vector<Real> & in,
                             Element /*global_element_id*/) override {
      auto nrows = spatial_dimension;
      auto nb_data = in.size() / (nrows * nrows);

      Vector<Real> von_mises_stress(nb_data);
      Matrix<Real> deviatoric_stress(3, 3);

      for (Int d = 0; d < nb_data; ++d) {
        MatrixProxy<const Real> cauchy_stress(in.data() + d * nrows * nrows,
                                              nrows, nrows);
        von_mises_stress(d) = Material::stressToVonMises(cauchy_stress);
      }

      return von_mises_stress;
    }

    Int getDim() override { return 1; }
    Int getNbComponent(Int /*old_nb_comp*/) override { return this->getDim(); }
  };

  /* ------------------------------------------------------------------------ */

} // namespace dumpers
} // namespace akantu

#endif /* AKANTU_DUMPER_MATERIAL_PADDERS_HH_ */
