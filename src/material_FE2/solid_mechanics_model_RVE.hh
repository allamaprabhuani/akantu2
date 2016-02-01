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
#include "solid_mechanics_model.hh"
#include "aka_grid_dynamic.hh"
#include <unordered_set>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__


class SolidMechanicsModelRVE : public SolidMechanicsModel {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  SolidMechanicsModelRVE(Mesh & mesh, bool use_RVE_mat_selector = true,
			 UInt nb_gel_pockets = 400,
			 UInt spatial_dimension = _all_dimensions,
			 const ID & id = "solid_mechanics_model",
			 const MemoryID & memory_id = 0);

  virtual ~SolidMechanicsModelRVE();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize completely the model
  virtual void initFull(const ModelOptions & options = SolidMechanicsModelOptions(_static, true));

  /// initialize the materials
  virtual void initMaterials();

  /// apply boundary contions based on macroscopic deformation gradient
  virtual void applyBoundaryConditions(const Matrix<Real> & displacement_gradient);

  /// advance the reactions -> grow gel and apply homogenized properties
  void advanceASR(const Matrix<Real> & prestrain);

  /// compute average stress or strain in the model
  Real averageTensorField(UInt row_index, UInt col_index, const ID & field_type);

  /// compute effective stiffness of the RVE
  void homogenizeStiffness(Matrix<Real> & C_macro);

  /// compute average eigenstrain
  void homogenizeEigenGradU(Matrix<Real> & eigen_gradu_macro);

  /// initialize the solver and the jacobian_matrix (called by initImplicit)
  virtual void initSolver(SolverOptions & options = _solver_no_options);

  /// allocate all vectors
  virtual void initArrays();

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */

  inline virtual void unpackData(CommunicationBuffer & buffer,
				 const UInt            index,
				 SynchronizationTag    tag);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(CornerNodes, corner_nodes, const Array<UInt> &);
  AKANTU_GET_MACRO(Volume, volume, Real);

private:

  /// find the corner nodes
  void findCornerNodes();

  /// perform virtual testing
  void performVirtualTesting(const Matrix<Real> & H, Matrix<Real> & eff_stresses, Matrix<Real> & eff_strains, const UInt test_no);

  /* ------------------------------------------------------------------------ */
  /* Members                                                                    */
  /* ------------------------------------------------------------------------ */

  /// volume of the RVE
  Real volume;

  /// corner nodes 1, 2, 3, 4 (see Leonardo's thesis, page 98)
  Array<UInt> corner_nodes;

  /// bottom nodes
  std::unordered_set<UInt> bottom_nodes;

  /// left nodes
  std::unordered_set<UInt> left_nodes;

  /// standard mat selector or user one
  bool use_RVE_mat_selector;

  StaticCommunicator * static_communicator_dummy;

  /// the number of gel pockets inside the RVE
  UInt nb_gel_pockets;
};


inline void SolidMechanicsModelRVE::unpackData(CommunicationBuffer & buffer,
					       const UInt index,
					       SynchronizationTag tag) {
  SolidMechanicsModel::unpackData(buffer, index, tag);

  if (tag == _gst_smm_uv) {
    Array<Real>::vector_iterator disp_it
      = displacement->begin(spatial_dimension);

    Vector<Real> current_disp(disp_it[index]);

    // if node is at the bottom, u_bottom = u_top +u_2 -u_3
    if ( bottom_nodes.count(index) ) {
      current_disp += Vector<Real>(disp_it[corner_nodes(1)]);
      current_disp -= Vector<Real>(disp_it[corner_nodes(2)]);
    }
    // if node is at the left, u_left = u_right +u_4 -u_3
    else if ( left_nodes.count(index) ) {
      current_disp += Vector<Real>(disp_it[corner_nodes(3)]);
      current_disp -= Vector<Real>(disp_it[corner_nodes(2)]);
    }
  }
}

/* -------------------------------------------------------------------------- */
/* ASR material selector                                                      */
/* -------------------------------------------------------------------------- */
class GelMaterialSelector : public MeshDataMaterialSelector<std::string>  {
public:
  GelMaterialSelector(SolidMechanicsModel & model,
		      const Real box_size,
		      const std::string & gel_material,
		      const UInt nb_gel_pockets,
		      Real tolerance = 0.) :
    MeshDataMaterialSelector<std::string>("physical_names", model),
    model(model),
    gel_material(gel_material),  
    nb_gel_pockets(nb_gel_pockets),
    nb_placed_gel_pockets(0),
    box_size(box_size) {
 
    Mesh & mesh = this->model.getMesh();
    UInt spatial_dimension = model.getSpatialDimension();

    const Vector<Real> & lower_bounds = mesh.getLowerBounds();
    const Vector<Real> & upper_bounds = mesh.getUpperBounds();

    Vector<Real> gcenter(spatial_dimension);
    for (UInt i = 0; i < spatial_dimension; ++i) {
      gcenter[i] = (upper_bounds[i] + lower_bounds[i]) / 2.;
    }

    Real grid_box_size = upper_bounds[0] - lower_bounds[0];
    Real grid_spacing = grid_box_size/7;
    Vector<Real> gspacing(spatial_dimension, grid_spacing);

    SpatialGrid<Element> grid(spatial_dimension, gspacing, gcenter);


    ElementType type = _triangle_3;
    GhostType ghost_type = _not_ghost;
    UInt nb_element = mesh.getNbElement(type, ghost_type);
    Element el;
    el.type = type;
    el.ghost_type = ghost_type;
    Array<Real> barycenter(0,2);
    barycenter.resize(nb_element);
    Array<Real>::vector_iterator bary_it = barycenter.begin(spatial_dimension);
    for (UInt elem = 0; elem < nb_element; ++bary_it, ++elem) {
      mesh.getBarycenter(elem, type, bary_it->storage(), ghost_type);
      el.element = elem;
      grid.insert(el, *bary_it);
    }

    /// generate the gel pockets
    srand(0.);
    Vector<Real> center(model.getSpatialDimension());
    for (UInt i = 0; i < this->nb_gel_pockets; ++i) {
      center.clear();
      center(0) = -box_size/2. + (box_size) * ((Real) rand() / (RAND_MAX));
      center(1) = -box_size/2. + (box_size) * ((Real) rand() / (RAND_MAX));
      Real min_dist = box_size;
      el.element = 0;
      bary_it = barycenter.begin(spatial_dimension);
      /// find cell in which current bary center lies
      SpatialGrid<Element>::CellID cell_id = grid.getCellID(center);
      SpatialGrid<Element>::Cell::const_iterator first_el =
	grid.beginCell(cell_id);
      SpatialGrid<Element>::Cell::const_iterator last_el =
	grid.endCell(cell_id);
      /// loop over all the elements in that cell
      for (;first_el != last_el; ++first_el){
	const Element & elem = *first_el;
	Vector<Real> bary = bary_it[elem.element];
	if (center.distance(bary) <= min_dist) {
	  min_dist = center.distance(bary);
	  el.element = elem.element;
	}
      }
      gel_pockets.push_back(el);
    }
  }

  UInt operator()(const Element & elem) {
    UInt temp_index = MeshDataMaterialSelector<std::string>::operator()(elem);
    if (temp_index == 1)
      return temp_index;
    std::vector<Element>::const_iterator iit = gel_pockets.begin();
    std::vector<Element>::const_iterator eit = gel_pockets.end();
    if(std::find(iit, eit, elem) != eit) {
      nb_placed_gel_pockets += 1;
      std::cout << nb_placed_gel_pockets << " gelpockets placed" << std::endl;
      return model.getMaterialIndex(gel_material);;
    }
    return 0;
  }


protected:
  SolidMechanicsModel & model;
  std::string gel_material;
  std::vector<Element> gel_pockets;
  UInt nb_gel_pockets;
  UInt nb_placed_gel_pockets;
  Real box_size;
};


__END_AKANTU__


///#include "material_selector_tmpl.hh"

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_RVE_HH__ */

