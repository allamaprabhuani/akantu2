/**
 * @file   asr_tools.hh
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Tue Jan 16 10:26:53 2014
 *
 * @brief  tools for the analysis of ASR samples
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
#include "aka_array.hh"
//#include "material_selector_tmpl.hh"
#include "mesh.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <unordered_set>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ASR_TOOLS_HH__
#define __AKANTU_ASR_TOOLS_HH__

namespace akantu {
class NodeGroup;
class SolidMechanicsModel;
} // namespace akantu

namespace akantu {

class ASRTools {
public:
  ASRTools(SolidMechanicsModel & model);

  virtual ~ASRTools() = default;

  /// This function is used in case of uni-axial boundary conditions:
  /// It detects the nodes, on which a traction needs to be applied and stores
  /// them in a node group
  void fillNodeGroup(NodeGroup & node_group, bool multi_axial = false);

  /// Apply boundary conditions on ASR samples to imitate lab testing conditions
  template <UInt dim>
  void applyBoundaryConditions(const bool free_expansion,
                               const Matrix<Real> & traction,
                               const ID & element_group, bool multi_axial);

  /// This function computes the average displacement along one side of the
  /// sample,
  /// caused by the expansion of gel pockets
  Real computeAverageDisplacement(SpatialDirection direction);

  /// This function computes a compoment of average volumetric strain tensor,
  /// i.e. eps_macro_xx or eps_macro_yy or eps_macro_zz
  Real computeVolumetricExpansion(SpatialDirection direction);

  /// This function computes the amount of damage in one material phase
  Real computeDamagedVolume(const ID & mat_name);

  /// This function is used to compute the average stiffness by performing a
  /// virtual tension test
  Real performTensionTest(SpatialDirection direction);

  /// This function calls all the functions to compute average properties and
  /// does the tension test
  void computeAveragePropertiesAndResidual(std::ofstream & file_output,
                                           Real time);

  /// perform tension tests and integrate the internal force on the upper
  /// surface
  void computeStiffnessReduction(std::ofstream & file_output, Real time);

  /// just the average properties, NO tension test
  void computeAverageProperties(std::ofstream & file_output);

  /// Same as the last one but writes a csv with first column having time in
  /// days
  void computeAverageProperties(std::ofstream & file_output, Real time);

  /// Save the state of the simulation for subsequent restart
  void saveState(UInt current_step);

  /// Load to previous state of the simulation for restart
  bool loadState(UInt & current_step);

  /// compute the volume occupied by a given material
  Real computePhaseVolume(const ID & mat_name);

  /// apply eigenstrain in cracks, that are filled with gel
  void applyEigenGradUinCracks(const Matrix<Real> prescribed_eigen_grad_u,
                               const ID & material_name);

  /// apply eigenstrain in the cracks, that advanced in the last step
  void applyEigenGradUinCracks(const Matrix<Real> prescribed_eigen_grad_u,
                               const ElementTypeMapUInt & critical_elements,
                               bool to_add = false);

  /// fill fully cracked elements with gel
  void fillCracks(ElementTypeMapReal & saved_damage);

  /// drain the gel in fully cracked elements
  void drainCracks(const ElementTypeMapReal & saved_damage);

  template <UInt dim> Real computeSmallestElementSize();

  // /// apply homogeneous temperature on Solid Mechanics Model
  // void applyTemperatureFieldToSolidmechanicsModel(const Real & temperature);

  /// compute increase in gel strain within 1 timestep
  Real computeDeltaGelStrainThermal(const Real delta_time, const Real k,
                             const Real activ_energy, const Real R,
                             const Real T, Real & amount_reactive_particles,
                             const Real saturation_const);

  /// compute linear increase in gel strain
  Real computeDeltaGelStrainLinear(const Real delta_time, const Real k);

  /* ------------------------------------------------------------------------ */
  /// RVE part

  /// apply boundary contions based on macroscopic deformation gradient
  virtual void
  applyBoundaryConditionsRve(const Matrix<Real> & displacement_gradient);

  /// apply homogeneous temperature field from the macroscale level to the RVEs
  virtual void applyHomogeneousTemperature(const Real & temperature);

  /// remove temperature from RVE on the end of ASR advancement
  virtual void removeTemperature();

  /// compute average stress or strain in the model
  Real averageTensorField(UInt row_index, UInt col_index,
                          const ID & field_type);

  /// compute effective stiffness of the RVE
  void homogenizeStiffness(Matrix<Real> & C_macro, bool first_time = false);

  /// compute average eigenstrain
  void homogenizeEigenGradU(Matrix<Real> & eigen_gradu_macro);

  /// compute damage volume in different phases
  void computeDamageRatio(Real & damage);

  /// dump the RVE
  void dumpRve();

private:
  /// find the corner nodes
  void findCornerNodes();

  /// perform virtual testing
  void performVirtualTesting(const Matrix<Real> & H,
                             Matrix<Real> & eff_stresses,
                             Matrix<Real> & eff_strains, const UInt test_no);

  // void fillCracks(ElementTypeMapReal & saved_damage);
  // void drainCracks(const ElementTypeMapReal & saved_damage);

  /* ------------------------------------------------------------------------ */
  /* Members */
  /* ------------------------------------------------------------------------ */

private:
  SolidMechanicsModel & model;

  // /// 2D hardcoded - no 3D support currently
  // using voigt_h = VoigtHelper<2>;

protected:
  /// volume of the RVE
  Real volume;

  /// corner nodes 1, 2, 3, 4 (see Leonardo's thesis, page 98)
  Array<UInt> corner_nodes;

  /// bottom nodes
  std::unordered_set<UInt> bottom_nodes;

  /// left nodes
  std::unordered_set<UInt> left_nodes;

  /// lower limit for stresses
  Array<Real> stress_limit;

  /// dump counter
  UInt nb_dumps;
};

/* -------------------------------------------------------------------------- */
/* ASR material selector                                                      */
/* -------------------------------------------------------------------------- */
class GelMaterialSelector : public MeshDataMaterialSelector<std::string> {
public:
  GelMaterialSelector(SolidMechanicsModel & model, const Real box_size,
                      const std::string & gel_material,
                      const UInt nb_gel_pockets,
                      std::string paste_material = "paste",
                      Real /*tolerance*/ = 0.)
      : MeshDataMaterialSelector<std::string>("physical_names", model),
        model(model), gel_material(gel_material),
        nb_gel_pockets(nb_gel_pockets), nb_placed_gel_pockets(0),
        box_size(box_size), paste_material(paste_material) {}

  void initGelPocket() {
    paste_material_id = model.getMaterialIndex(paste_material);

    Mesh & mesh = this->model.getMesh();
    UInt dim = model.getSpatialDimension();
    //    Element el{_triangle_3, 0, _not_ghost};
    for (auto el_type :
         model.getMaterial("aggregate").getElementFilter().elementTypes(dim)) {

      const auto & filter =
          model.getMaterial("aggregate").getElementFilter()(el_type);
      if (!filter.size() == 0)
        AKANTU_EXCEPTION("Check the element type for aggregate material");

      Element el{el_type, 0, _not_ghost};
      UInt nb_element = mesh.getNbElement(el.type, el.ghost_type);
      Array<Real> barycenter(nb_element, dim);

      for (auto && data : enumerate(make_view(barycenter, dim))) {
        el.element = std::get<0>(data);
        auto & bary = std::get<1>(data);
        mesh.getBarycenter(el, bary);
      }

      /// generate the gel pockets
      srand(0.);
      Vector<Real> center(dim);
      UInt placed_gel_pockets = 0;
      std::set<int> checked_baries;
      while (placed_gel_pockets != nb_gel_pockets) {
        /// get a random bary center
        UInt bary_id = rand() % nb_element;
        if (checked_baries.find(bary_id) != checked_baries.end())
          continue;
        checked_baries.insert(bary_id);
        el.element = bary_id;
        if (MeshDataMaterialSelector<std::string>::operator()(el) ==
            paste_material_id)
          continue; /// element belongs to paste
        gel_pockets.push_back(el);
        placed_gel_pockets += 1;
      }
    }
    is_gel_initialized = true;
  }

  UInt operator()(const Element & elem) {
    if (not is_gel_initialized)
      initGelPocket();

    UInt temp_index = MeshDataMaterialSelector<std::string>::operator()(elem);
    if (temp_index == paste_material_id)
      return temp_index;
    auto iit = gel_pockets.begin();
    auto eit = gel_pockets.end();
    if (std::find(iit, eit, elem) != eit) {
      nb_placed_gel_pockets += 1;
      std::cout << nb_placed_gel_pockets << " gelpockets placed" << std::endl;
      return model.getMaterialIndex(gel_material);
      ;
    }
    return temp_index;
  }

protected:
  SolidMechanicsModel & model;
  std::string gel_material;
  std::vector<Element> gel_pockets;
  UInt nb_gel_pockets;
  UInt nb_placed_gel_pockets;
  Real box_size;
  std::string paste_material{"paste"};
  UInt paste_material_id{1};
  bool is_gel_initialized{false};
};

} // namespace akantu

#endif /* __AKANTU_ASR_TOOLS_HH__ */
