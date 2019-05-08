/**
 * @file   model_options.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 04 2017
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  A Documented file.
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_named_argument.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MODEL_OPTIONS_HH__
#define __AKANTU_MODEL_OPTIONS_HH__

namespace akantu {

namespace {
  DECLARE_NAMED_ARGUMENT(analysis_method);
}

struct ModelOptions {
  explicit ModelOptions(AnalysisMethod analysis_method = _static)
      : analysis_method(analysis_method) {}

  template <typename... pack>
  ModelOptions(use_named_args_t, pack &&... _pack)
      : ModelOptions(OPTIONAL_NAMED_ARG(analysis_method, _static)) {}

  virtual ~ModelOptions() = default;

  AnalysisMethod analysis_method;
};

#ifdef AKANTU_SOLID_MECHANICS
/* -------------------------------------------------------------------------- */
struct SolidMechanicsModelOptions : public ModelOptions {
  explicit SolidMechanicsModelOptions(
      AnalysisMethod analysis_method = _explicit_lumped_mass)
      : ModelOptions(analysis_method) {}

  template <typename... pack>
  SolidMechanicsModelOptions(use_named_args_t, pack &&... _pack)
      : SolidMechanicsModelOptions(
            OPTIONAL_NAMED_ARG(analysis_method, _explicit_lumped_mass)) {}
};
#endif

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_COHESIVE_ELEMENT
namespace {
  DECLARE_NAMED_ARGUMENT(is_extrinsic);
}
/* -------------------------------------------------------------------------- */
struct SolidMechanicsModelCohesiveOptions : public SolidMechanicsModelOptions {
  SolidMechanicsModelCohesiveOptions(
      AnalysisMethod analysis_method = _explicit_lumped_mass,
      bool extrinsic = false)
      : SolidMechanicsModelOptions(analysis_method), is_extrinsic(extrinsic) {}

  template <typename... pack>
  SolidMechanicsModelCohesiveOptions(use_named_args_t, pack &&... _pack)
      : SolidMechanicsModelCohesiveOptions(
            OPTIONAL_NAMED_ARG(analysis_method, _explicit_lumped_mass),
            OPTIONAL_NAMED_ARG(is_extrinsic, false)) {}

  bool is_extrinsic{false};
};
#endif

#ifdef AKANTU_HEAT_TRANSFER
/* -------------------------------------------------------------------------- */
struct HeatTransferModelOptions : public ModelOptions {
  explicit HeatTransferModelOptions(
      AnalysisMethod analysis_method = _explicit_lumped_mass)
      : ModelOptions(analysis_method) {}

  template <typename... pack>
  HeatTransferModelOptions(use_named_args_t, pack &&... _pack)
      : HeatTransferModelOptions(
            OPTIONAL_NAMED_ARG(analysis_method, _explicit_lumped_mass)) {}
};
#endif

#ifdef AKANTU_EMBEDDED

namespace {
  DECLARE_NAMED_ARGUMENT(init_intersections);
}
/* -------------------------------------------------------------------------- */
struct EmbeddedInterfaceModelOptions : SolidMechanicsModelOptions {
  /**
   * @brief Constructor for EmbeddedInterfaceModelOptions
   * @param analysis_method see SolidMechanicsModelOptions
   * @param init_intersections compute intersections
   */
  EmbeddedInterfaceModelOptions(
      AnalysisMethod analysis_method = _explicit_lumped_mass,
      bool init_intersections = true)
      : SolidMechanicsModelOptions(analysis_method),
        has_intersections(init_intersections) {}

  template <typename... pack>
  EmbeddedInterfaceModelOptions(use_named_args_t, pack &&... _pack)
      : EmbeddedInterfaceModelOptions(
            OPTIONAL_NAMED_ARG(analysis_method, _explicit_lumped_mass),
            OPTIONAL_NAMED_ARG(init_intersections, true)) {}

  /// Should consider reinforcements
  bool has_intersections;
};
#endif


#ifdef AKANTU_CONTACT_MECHANICS
namespace {
  DECLARE_NAMED_ARGUMENT(is_explicit);
}
/* -------------------------------------------------------------------------- */
struct ContactMechanicsModelOptions : public ModelOptions {
  explicit ContactMechanicsModelOptions(
      AnalysisMethod analysis_method = _explicit_contact,
      bool explicit = true)
    : ModelOptions(analysis_method), is_explicit(explicit) {}

  template <typename... pack>
  ContactMechanicsModelOptions(use_named_args_t, pack &&... _pack)
      : ContactMechanicsModelOptions(
	      OPTIONAL_NAMED_ARG(analysis_method, _explicit_contact),
	      OPTIONAL_NAMED_ARG(is_explicit, true)) {}

    bool is_explicit{true};
};
#endif

#ifdef AKANTU_MODEL_COUPLERS
namespace {
  DECLARE_NAMED_ARGUMENT(is_explicit);
}
/* -------------------------------------------------------------------------- */
struct CouplerSolidContactOptions : public ModelOptions {
  explicit CouplerSolidContactOptions(
      AnalysisMethod analysis_method = _explicit_contact,
      bool explicit = true)
    : ModelOptions(analysis_method), is_explicit(explicit) {}

  template <typename... pack>
  CouplerSolidContactOptions(use_named_args_t, pack &&... _pack)
      : CouplerSolidContactOptions(
	      OPTIONAL_NAMED_ARG(analysis_method, _explicit_contact),
	      OPTIONAL_NAMED_ARG(is_explicit, true)) {}

    bool is_explicit{true};
};  
#endif
  
} // akantu

#endif /* __AKANTU_MODEL_OPTIONS_HH__ */
