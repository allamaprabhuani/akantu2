/**
 * @file   material_reinforcement.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Mar 13 2015
 * @date last modification: Tue Nov 24 2015
 *
 * @brief  Reinforcement material
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

#ifndef __AKANTU_MATERIAL_REINFORCEMENT_HH__
#define __AKANTU_MATERIAL_REINFORCEMENT_HH__

#include "aka_common.hh"

#include "material.hh"
#include "embedded_interface_model.hh"
#include "embedded_internal_field.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/**
 * @brief Material used to represent embedded reinforcements
 *
 * This class is used for computing the reinforcement stiffness matrix
 * along with the reinforcement residual. Room is made for constitutive law,
 * but actual use of contitutive laws is made in MaterialReinforcementTemplate.
 *
 * Be careful with the dimensions in this class :
 *  -  this->spatial_dimension is always 1
 *  -  the template parameter dim is the dimension of the problem
 */
template<UInt dim>
class MaterialReinforcement : virtual public Material {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /// Constructor
  MaterialReinforcement(SolidMechanicsModel & model,
                        UInt spatial_dimension,
                        const Mesh & mesh,
                        FEEngine & fe_engine,
                        const ID & id = "");

  /// Destructor
  virtual ~MaterialReinforcement();

protected:
  void initialize(SolidMechanicsModel & a_model);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// Init the material
  virtual void initMaterial();

  /// Init the background shape derivatives
  void initBackgroundShapeDerivatives();

  /// Init the cosine matrices
  void initDirectingCosines();

  /// Assemble stiffness matrix
  virtual void assembleStiffnessMatrix(GhostType ghost_type);

  /// Update the residual
  virtual void updateResidual(GhostType ghost_type = _not_ghost);

  /// Assembled the residual
  virtual void assembleResidual(GhostType ghost_type);

  /// Compute all the stresses !
  virtual void computeAllStresses(GhostType ghost_type);

  /// Compute the stiffness parameter for elements of a type
  virtual void computeTangentModuli(const ElementType & type,
                                    Array<Real> & tangent,
                                    GhostType ghost_type) = 0;

  /// Compute energy
  virtual Real getEnergy(std::string id);

  // virtual ElementTypeMap<UInt> getInternalDataPerElem(const ID & field_name,
  //                                                     const ElementKind & kind,
  //                                                     const ID & fe_engine_id) const;

  // /// Reimplementation of Material's function to accomodate for interface mesh
  // virtual void flattenInternal(const std::string & field_id,
  //                              ElementTypeMapArray<Real> & internal_flat,
  //                              const GhostType ghost_type = _not_ghost,
  //                              ElementKind element_kind = _ek_not_defined) const;

  /* ------------------------------------------------------------------------ */
  /* Protected methods                                                        */
  /* ------------------------------------------------------------------------ */
protected:
  /// Allocate the background shape derivatives
  void allocBackgroundShapeDerivatives();

  /// Compute the directing cosines matrix for one element type
  void computeDirectingCosines(const ElementType & type, GhostType ghost_type);

  /// Compute the directing cosines matrix on quadrature points.
  inline void computeDirectingCosinesOnQuad(const Matrix<Real> & nodes,
                                            Matrix<Real> & cosines);

  /// Assemble the stiffness matrix for an element type (typically _segment_2)
  void assembleStiffnessMatrix(const ElementType & type, GhostType ghost_type);

  /// Assemble the stiffness matrix for background & interface types
  void assembleStiffnessMatrixInterface(const ElementType & interface_type,
                                        const ElementType & background_type,
                                        GhostType ghost_type);

  /// Compute the background shape derivatives for a type
  void computeBackgroundShapeDerivatives(const ElementType & type, GhostType ghost_type);

  /// Filter elements crossed by interface of a type
  void filterInterfaceBackgroundElements(Array<UInt> & filter,
                                         const ElementType & type,
                                         const ElementType & interface_type,
                                         GhostType ghost_type);

  /// Assemble the residual of one type of element (typically _segment_2)
  void assembleResidual(const ElementType & type, GhostType ghost_type);

  /// Assemble the residual for a pair of elements
  void assembleResidualInterface(const ElementType & interface_type,
                                 const ElementType & background_type,
                                 GhostType ghost_type);

  // TODO figure out why voigt size is 4 in 2D
  inline void stressTensorToVoigtVector(const Matrix<Real> & tensor, Vector<Real> & vector);
  inline void strainTensorToVoigtVector(const Matrix<Real> & tensor, Vector<Real> & vector);

  /// Compute gradu on the interface quadrature points
  virtual void computeGradU(const ElementType & type, GhostType ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Embedded model
  EmbeddedInterfaceModel * model;

  /// Stress in the reinforcement
  InternalField<Real> stress_embedded;

  /// Gradu of concrete on reinforcement
  InternalField<Real> gradu_embedded;

  /// C matrix on quad
  InternalField<Real> directing_cosines;

  /// Prestress on quad
  InternalField<Real> pre_stress;

  /// Cross-sectional area
  Real area;

  /// Background mesh shape derivatives
  ElementTypeMap< ElementTypeMapArray<Real> * > shape_derivatives;

};

#include "material_reinforcement_inline_impl.cc"

__END_AKANTU__

#endif // __AKANTU_MATERIAL_REINFORCEMENT_HH__
