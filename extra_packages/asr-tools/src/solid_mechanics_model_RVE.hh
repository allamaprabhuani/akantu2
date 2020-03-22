/**
 * @file   solid_mechanics_model_RVE.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Jan 13 14:54:18 2016
 *
 * @brief  SMM for RVE computations in FE2 simulations
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
#ifndef __AKANTU_SOLID_MECHANICS_MODEL_RVE_HH__
#define __AKANTU_SOLID_MECHANICS_MODEL_RVE_HH__

/* -------------------------------------------------------------------------- */
#include "aka_grid_dynamic.hh"
#include "asr_tools.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
//#include <unordered_set>
/* -------------------------------------------------------------------------- */

namespace akantu {

class SolidMechanicsModelRVE : public SolidMechanicsModel, public ASRTools {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  SolidMechanicsModelRVE(Mesh & mesh, bool use_RVE_mat_selector = true,
                         UInt nb_gel_pockets = 400, UInt dim = _all_dimensions,
                         const ID & id = "solid_mechanics_model",
                         const MemoryID & memory_id = 0);

  virtual ~SolidMechanicsModelRVE();

  typedef VoigtHelper<2> voigt_h;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  void initFullImpl(const ModelOptions & option) override;

  /// initialize the materials
  void initMaterials() override;

public:
  /// advance the reactions -> grow gel and apply homogenized properties
  void advanceASR(const Matrix<Real> & prestrain);

  /// correct the rigid boundary movement and assemble internal forces
  void assembleInternalForces() override;
  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<UInt> & index,
                         const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(CornerNodes, corner_nodes, const Array<UInt> &);
  AKANTU_GET_MACRO(Volume, volume, Real);
  bool hasStiffnessChanged() { return this->stiffness_changed; };

private:
  /* ------------------------------------------------------------------------ */
  /* Members */
  /* ------------------------------------------------------------------------ */
  /// standard mat selector or user one
  bool use_RVE_mat_selector;

  /// the number of gel pockets inside the RVE
  UInt nb_gel_pockets;

  /// the number of gel pockets inside the RVE
  bool stiffness_changed;
};

inline void SolidMechanicsModelRVE::unpackData(CommunicationBuffer & buffer,
                                               const Array<UInt> & index,
                                               const SynchronizationTag & tag) {
  SolidMechanicsModel::unpackData(buffer, index, tag);

  //  if (tag == _gst_smm_uv) {
  //    auto disp_it = displacement->begin(spatial_dimension);
  //
  //    for (auto node : index) {
  //      Vector<Real> current_disp(disp_it[node]);
  //
  //      // if node is at the bottom, u_bottom = u_top +u_2 -u_3
  //      if (bottom_nodes.count(node)) {
  //        current_disp += Vector<Real>(disp_it[corner_nodes(1)]);
  //        current_disp -= Vector<Real>(disp_it[corner_nodes(2)]);
  //      }
  //      // if node is at the left, u_left = u_right +u_4 -u_3
  //      else if (left_nodes.count(node)) {
  //        current_disp += Vector<Real>(disp_it[corner_nodes(3)]);
  //        current_disp -= Vector<Real>(disp_it[corner_nodes(2)]);
  //      }
  //    }
  //  }
}

// /* --------------------------------------------------------------------------
// */
// /* ASR material selector */
// /* --------------------------------------------------------------------------
// */ class GelMaterialSelector : public MeshDataMaterialSelector<std::string> {
// public:
//   GelMaterialSelector(SolidMechanicsModel & model, const Real box_size,
//                       const std::string & gel_material,
//                       const UInt nb_gel_pockets,
//                       std::string paste_material = "paste",
//                       Real /*tolerance*/ = 0.)
//       : MeshDataMaterialSelector<std::string>("physical_names", model),
//         model(model), gel_material(gel_material),
//         nb_gel_pockets(nb_gel_pockets), nb_placed_gel_pockets(0),
//         box_size(box_size), paste_material(paste_material) {}

//   void initGelPocket() {
//     paste_material_id = model.getMaterialIndex(paste_material);

//     Mesh & mesh = this->model.getMesh();
//     UInt dim = model.getSpatialDimension();
//     //    Element el{_triangle_3, 0, _not_ghost};
//     for (auto el_type :
//          model.getMaterial("aggregate").getElementFilter().elementTypes(dim))
//          {

//       const auto & filter =
//       model.getMaterial("aggregate").getElementFilter()(el_type); if
//       (!filter.size() == 0)
//         AKANTU_EXCEPTION("Check the element type for aggregate material");

//       Element el{el_type, 0, _not_ghost};
//       UInt nb_element = mesh.getNbElement(el.type, el.ghost_type);
//       Array<Real> barycenter(nb_element, dim);

//       for (auto && data : enumerate(make_view(barycenter, dim))) {
//         el.element = std::get<0>(data);
//         auto & bary = std::get<1>(data);
//         mesh.getBarycenter(el, bary);
//       }

//       /// generate the gel pockets
//       srand(0.);
//       Vector<Real> center(dim);
//       UInt placed_gel_pockets = 0;
//       std::set<int> checked_baries;
//       while (placed_gel_pockets != nb_gel_pockets) {
//         /// get a random bary center
//         UInt bary_id = rand() % nb_element;
//         if (checked_baries.find(bary_id) != checked_baries.end())
//           continue;
//         checked_baries.insert(bary_id);
//         el.element = bary_id;
//         if (MeshDataMaterialSelector<std::string>::operator()(el) ==
//             paste_material_id)
//           continue; /// element belongs to paste
//         gel_pockets.push_back(el);
//         placed_gel_pockets += 1;
//       }
//     }
//     is_gel_initialized = true;
//   }

//   UInt operator()(const Element & elem) {
//     if (not is_gel_initialized)
//       initGelPocket();

//     UInt temp_index =
//     MeshDataMaterialSelector<std::string>::operator()(elem); if (temp_index
//     == paste_material_id)
//       return temp_index;
//     auto iit = gel_pockets.begin();
//     auto eit = gel_pockets.end();
//     if (std::find(iit, eit, elem) != eit) {
//       nb_placed_gel_pockets += 1;
//       std::cout << nb_placed_gel_pockets << " gelpockets placed" <<
//       std::endl; return model.getMaterialIndex(gel_material);
//       ;
//     }
//     return temp_index;
//   }

// protected:
//   SolidMechanicsModel & model;
//   std::string gel_material;
//   std::vector<Element> gel_pockets;
//   UInt nb_gel_pockets;
//   UInt nb_placed_gel_pockets;
//   Real box_size;
//   std::string paste_material{"paste"};
//   UInt paste_material_id{1};
//   bool is_gel_initialized{false};
// }; // namespace akantu

} // namespace akantu

///#include "material_selector_tmpl.hh"

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_RVE_HH__ */
